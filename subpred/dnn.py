from tensorflow import keras
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import (
    RepeatedStratifiedKFold,
)


# Test result AT sugar amino: the three models perform very similarly, just use the simplest one (create_model)
# Tried different values for dropout (0.0,0.3,0.5,0.7), performance for non-0 is similar.
# Dropout 0.5 for smaller datasets, maybe try 0.3 for larger
# Tried different values for batch_size. Not really a big difference. 8 for smaller (less overfitting), 16 or 32 for larger (more speed)
# TODO put this code into subpred package, once TF is compatible with new GPU


def create_model(n_features):
    # Larger datasets: try lower dropout
    # Try starting at lower number of nodes
    model = keras.Sequential(
        [
            keras.layers.Input(shape=(n_features,)),
            keras.layers.Dense(512, activation="relu"),
            keras.layers.Dropout(0.5),
            keras.layers.Dense(256, activation="relu"),
            keras.layers.Dropout(0.5),
            keras.layers.Dense(128, activation="relu"),
            keras.layers.Dense(1, activation="sigmoid"),
        ]
    )

    model.compile(
        optimizer="adam",
        loss="binary_crossentropy",
        metrics=[
            keras.metrics.F1Score(average="macro", name="F1_macro"),
            keras.metrics.TruePositives(name="TP"),
            keras.metrics.TrueNegatives(name="TN"),
            keras.metrics.FalsePositives(name="FP"),
            keras.metrics.FalseNegatives(name="FN"),
        ],
    )
    return model


def create_model_dynamic_nodes(n_features):
    if n_features > 1024:
        layer_sizes = [1024, 512, 256]
        print("selecting large model")
    elif n_features > 512:  # includes embeddings with len 1024
        layer_sizes = [512, 256, 128]
        print("selecting medium model")
    else:
        layer_sizes = [256, 128, 64]
        print("selecting small model")
    model = keras.Sequential(
        [
            keras.layers.Input(shape=(n_features,)),
            keras.layers.Dense(layer_sizes[0], activation="relu"),
            keras.layers.Dropout(0.5),
            keras.layers.Dense(layer_sizes[1], activation="relu"),
            keras.layers.Dropout(0.5),
            keras.layers.Dense(layer_sizes[2], activation="relu"),
            keras.layers.Dense(1, activation="sigmoid"),
        ]
    )

    model.compile(
        optimizer="adam",
        loss="binary_crossentropy",
        metrics=[
            keras.metrics.F1Score(average="macro", name="F1_macro"),
            keras.metrics.TruePositives(name="TP"),
            keras.metrics.TrueNegatives(name="TN"),
            keras.metrics.FalsePositives(name="FP"),
            keras.metrics.FalseNegatives(name="FN"),
        ],
    )
    return model


def create_model_dynamic_layers(n_features):
    model = keras.Sequential()
    model.add(keras.layers.Input(shape=(n_features,)))
    for layer_size in [2048, 1024, 512, 256]:
        if n_features >= layer_size:
            model.add(keras.layers.Dense(layer_size, activation="relu"))
            model.add(keras.layers.Dropout(0.5))
    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dense(1, activation="sigmoid"))

    model.compile(
        optimizer="adam",
        loss="binary_crossentropy",
        metrics=[
            keras.metrics.F1Score(average="macro", name="F1_macro"),
            keras.metrics.TruePositives(name="TP"),
            keras.metrics.TrueNegatives(name="TN"),
            keras.metrics.FalsePositives(name="FP"),
            keras.metrics.FalseNegatives(name="FN"),
        ],
    )

    return model


from sklearn.utils.class_weight import compute_class_weight


def crossval_dnn(
    ml_dataset,
    model_func,
    scores_dict,
    splits=5,
    repeats=5,
    epochs=100,
    batch_size=8,
    verbose=False,
    calculate_class_weights=False,
):
    print(f"=== {ml_dataset.name} ===")
    preprocess = make_pipeline(VarianceThreshold(0.0), StandardScaler())

    X, y = ml_dataset.X, ml_dataset.y

    train_scores = list()
    test_scores = list()
    fold_count = 1
    for train_idx_outer, val_idx_outer in RepeatedStratifiedKFold(
        n_splits=splits, n_repeats=repeats, random_state=0
    ).split(X, y):
        if verbose:
            print(f"Fold {fold_count} out of {splits*repeats}")
        fold_count += 1

        X_train, X_test = X[train_idx_outer], X[val_idx_outer]
        y_train, y_test = y[train_idx_outer], y[val_idx_outer]

        X_train = preprocess.fit_transform(X_train, y_train)
        X_test = preprocess.transform(X_test)

        # important: create from scratch to reset weights
        model = model_func(X_train.shape[1])
        # TODO Early Stopping can be an option for larger datasets (less overfitting, faster training)
        # TODO class weights option
        class_weights_dict = None
        if calculate_class_weights:
            classes = np.sort(np.unique(y_train))
            class_weights = compute_class_weight(
                class_weight="balanced", classes=classes, y=y_train
            )
            class_weights_dict = dict(zip(classes, class_weights))
        training_history = model.fit(
            X_train,
            y_train.reshape(-1, 1),
            epochs=epochs,
            batch_size=batch_size,
            verbose="auto" if verbose else 0,
            class_weight=class_weights_dict,
        )
        y_prob = model.predict(X_test, verbose="auto" if verbose else 0)
        y_pred = (y_prob > 0.5).astype(int).flatten()
        # TODO log mis-classified samples: always the same ones?

        for score_name, score_func in scores_dict.items():
            test_scores.append(
                (ml_dataset.name, score_name, score_func(y_test, y_pred))
            )
            if verbose:
                print(score_name, score_func(y_test, y_pred))

        res = model.evaluate(
            X_test, y_test.reshape(-1, 1), verbose="auto" if verbose else 0
        )
        if verbose:
            print(res)
    df_scores = pd.DataFrame(test_scores, columns=["Feature", "Metric", "Value"])

    return df_scores

# Example dict:
# scoring_outer = {
#     "Balanced Accuracy": balanced_accuracy_score,
#     "F1 Macro": lambda y_test, y_pred: f1_score(
#         y_true=y_test, y_pred=y_pred, average="macro"
#     ),
#     "F1 Class 0": lambda y_test, y_pred: f1_score(
#         y_true=y_test, y_pred=y_pred, pos_label=0
#     ),
#     "F1 Class 1": lambda y_test, y_pred: f1_score(
#         y_true=y_test, y_pred=y_pred, pos_label=1
#     ),
# }