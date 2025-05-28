"""
@author: adenger
"""

from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.svm import LinearSVC, SVC
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import (
    GridSearchCV,
    cross_val_score,
    RepeatedStratifiedKFold,
    StratifiedKFold,
)

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from subpred.util import save_data, load_data
import matplotlib.pyplot as plt


# Tested: gives same results but without error messages
class DynamicSelectKBest(BaseEstimator, TransformerMixin):
    def __init__(self, score_func=f_classif, k=10):
        self.score_func = score_func
        self.k = k
        self.selector_ = None

    def fit(self, X, y=None):
        # Avoids error where variancethreshold removes features,
        # then selectkbest tries to select more features than are available.
        # Gridsearch sets value of k directly
        k = min(self.k, X.shape[1])
        self.selector_ = SelectKBest(score_func=self.score_func, k=k)
        self.selector_.fit(X, y)
        return self

    def transform(self, X):
        return self.selector_.transform(X)

    def get_support(self, indices=False):
        return self.selector_.get_support(indices=indices)


def nested_crossval_svm(
    test_name,
    X,
    y,
    sample_names,
    feature_names,
    outer_cv: int = 5,
    inner_cv: int = 5,
    repeats: int = 10,
    n_jobs_inner: int = 1,
    n_jobs_outer: int = -1,
    scoring_inner: str = "balanced_accuracy",
    scoring_outer: str = "balanced_accuracy",
):
    print(f"=== {test_name} ===")
    model = make_pipeline(
        VarianceThreshold(), StandardScaler(), DynamicSelectKBest(), SVC()
    )

    max_features = min(len(feature_names), 200)

    param_grid = {
        "dynamicselectkbest__k": list(range(1, max_features, 1)),
        "svc__class_weight": ["balanced"],
        "svc__C": [0.1, 1, 10],
        "svc__gamma": ["scale"],
    }
    # scale : 1 / (n_features * X.var()).
    # larger variance and more features leads to a smoother decision boundary,
    # where each sample has less influence

    gridsearch = GridSearchCV(
        estimator=model,
        param_grid=param_grid,
        scoring=scoring_inner,
        cv=StratifiedKFold(inner_cv),
        n_jobs=n_jobs_inner,
    )

    # Nested loop (sound results):
    # gridsearch.n_jobs = 1
    nested_crossval_results = cross_val_score(
        gridsearch,
        X,
        y,
        cv=RepeatedStratifiedKFold(
            n_splits=outer_cv, n_repeats=repeats, random_state=0
        ),
        scoring=scoring_outer,
        n_jobs=n_jobs_outer,
    )
    print(
        f"Nested crossvalidation: {nested_crossval_results.mean():.2f}+-{nested_crossval_results.std():.2f}"
    )
    return nested_crossval_results


def get_svm_results(
    ml_datasets: list,
    output_folder: str,
    test_name: str,
    recalculate: bool = True,
    **kwargs  # for nested_crossval
):
    # wrapper method for caching, and preparing long-form results for plot
    # performs nested crossval for each feature in ml_datasets
    if not recalculate and Path(output_folder + test_name + ".pickle").exists():
        df_results_long = load_data(test_name, folder_path=output_folder)
    else:

        results_rbf_svm = [
            (
                ml_dataset[0],  
                nested_crossval_svm(
                    *ml_dataset,
                    **kwargs
                ),
            )
            for ml_dataset in ml_datasets
        ]

        results_long = list()
        for feature_name, test_results in results_rbf_svm:
            for test_result in test_results:
                results_long.append((feature_name, test_result))

        if "scoring_outer" in kwargs:
            score_name = kwargs["scoring_outer"].replace("_"," ").title()
        else:
            score_name = "Score"

        df_results_long = pd.DataFrame.from_records(
            results_long, columns=["Feature Name", score_name]
        )

        save_data(df_results_long, test_name, folder_path=output_folder)
    return df_results_long


def plot_results_long(
    df_results_long: pd.DataFrame, output_folder_path: str, test_name: str
):
    plt.figure(figsize=(10, 5), dpi=500)
    score_name = df_results_long.columns[1]
    print(f"creating plot for score {score_name}")
    df_results_long_stats = pd.concat(
        [
            df_results_long.groupby("Feature Name")[score_name]
            .mean()  # TODO mean?
            .rename("mean_val"),
            df_results_long.groupby("Feature Name")[score_name]
            .std()
            .rename("std_val"),
        ],
        axis=1,
        verify_integrity=True
    )
    df_results_long_stats = df_results_long_stats.sort_values(
        by=["mean_val", "std_val"], ascending=[True, False]
    )
    sns.boxplot(
        df_results_long,
        x="Feature Name",
        y=score_name,
        order=df_results_long_stats.index,
    )
    df_results_long
    plt.xticks(rotation=90)
    plt.ylim((0, 1.05))
    plt.grid(True, alpha=0.5)
    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.savefig(output_folder_path + test_name, bbox_inches="tight", dpi=300)
    return df_results_long_stats


# TODO This whole cell as function


def find_outliers(X: np.array):
    outlier_detector = make_pipeline(
        StandardScaler(),
        PCA(n_components=0.95),
        IsolationForest(contamination="auto", random_state=0),
    )
    outliers = outlier_detector.fit_predict(X)  # -1 for outliers, 1 for inliers
    is_outlier = outliers == -1
    # print(is_outlier.sum(), "outliers found")
    return is_outlier


def outlier_check(dataset_full, ml_datasets: list, threshold: float = 0.8):

    outliers = np.array(
        [find_outliers(X) for feature_name, X, _, sample_names, _ in ml_datasets]
    )

    df_outliers = (
        pd.Series(index=ml_datasets[0][3], data=outliers.sum(axis=0))
        .sort_values(ascending=False)
        .to_frame(name="outlier_count")
        .join(dataset_full[0].protein_names)
        .join(dataset_full[1])
    )

    df_outliers = df_outliers[df_outliers.outlier_count >= len(ml_datasets) * threshold]
    # potential outliers to manually check TODO useful?
    return df_outliers
