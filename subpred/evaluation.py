"""
@author: adenger
"""

from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC, SVC
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import (
    GridSearchCV,
    cross_val_score,
    cross_validate,
    RepeatedStratifiedKFold,
    StratifiedKFold,
)

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from subpred.util import save_data, load_data
import matplotlib.pyplot as plt
import seaborn as sns
from subpred.mldataset import MLDataset


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
    ml_dataset: MLDataset,
    outer_cv: int = 5,
    inner_cv: int = 5,
    repeats: int = 1,
    n_jobs_inner: int = 1,
    n_jobs_outer: int = -1,
    scoring_inner: str = "balanced_accuracy",
    scoring_outer: dict = {"Balanced Accuracy": "balanced_accuracy"},
    svm_C1: bool = False,
    max_k: int = 200,
):
    print(f"=== {ml_dataset.name} ===")
    model = make_pipeline(
        VarianceThreshold(), StandardScaler(), DynamicSelectKBest(), SVC()
    )

    max_features = min(len(ml_dataset.feature_names), max_k)

    param_grid = {
        "dynamicselectkbest__k": list(range(1, max_features, 1)),
        "svc__class_weight": ["balanced"],
        "svc__C": [0.1, 1, 10],
        "svc__gamma": ["scale"],
    }
    if svm_C1:
        param_grid["svc__C"] = [1]
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
    nested_crossval_results = cross_validate(
        gridsearch,
        ml_dataset.X,
        ml_dataset.y,
        cv=RepeatedStratifiedKFold(
            n_splits=outer_cv, n_repeats=repeats, random_state=0
        ),
        scoring=scoring_outer,
        n_jobs=n_jobs_outer,
    )
    cv_results = [
        (metric_name, nested_crossval_results[f"test_{metric_name}"])
        for metric_name, _ in scoring_outer.items()
    ]
    # print(nested_crossval_results)
    # if isinstance(scoring_outer, str):
    #     cv_results = [(scoring_outer, nested_crossval_results["test_score"])]
    # else:  # list
    #     cv_results = [
    #         (metric_name, nested_crossval_results[f"test_{metric_name}"])
    #         for metric_name in scoring_outer
    #     ]

    for score_name, scores in cv_results:
        print(f"{score_name}: {scores.mean():.2f}+-{scores.std():.2f}")
    return cv_results


def get_svm_results(
    ml_datasets: list,
    output_folder: str,
    test_name: str,
    recalculate: bool = True,
    **kwargs,  # for nested_crossval
):
    # wrapper method for caching, and preparing long-form results for plot
    # performs nested crossval for each feature in ml_datasets
    if not recalculate and Path(output_folder + test_name + ".pickle").exists():
        df_results_long = load_data(test_name, folder_path=output_folder)
    else:

        results_rbf_svm = [
            (
                ml_dataset.name,
                nested_crossval_svm(ml_dataset, **kwargs),
            )
            for ml_dataset in ml_datasets
        ]

        results_long = list()
        for feature_name, test_results in results_rbf_svm:
            for metric_name, metric_scores in test_results:
                for metric_score in metric_scores:
                    results_long.append((feature_name, metric_name, metric_score))

        # if "scoring_outer" in kwargs:
        #     score_name = kwargs["scoring_outer"].replace("_"," ").title()
        # else:
        #     score_name = "Score"

        df_results_long = pd.DataFrame.from_records(
            results_long, columns=["Feature", "Metric", "Value"]
        )

        save_data(df_results_long, test_name, folder_path=output_folder)
    return df_results_long


def summarize_results_long(df_results_long: pd.DataFrame):
    return pd.concat(
        [
            df_results_long.groupby(["Feature", "Metric"])
            .mean()
            .rename(columns={"Value": "Mean"}),
            df_results_long.groupby(["Feature", "Metric"])
            .median()
            .rename(columns={"Value": "Median"}),
            df_results_long.groupby(["Feature", "Metric"])
            .std()
            .rename(columns={"Value": "Sdev"}),
        ],
        axis=1,
    )


def plot_results_long(
    df_results_long: pd.DataFrame,
    output_folder_path: str,
    test_name: str,
    plot_order: list = None,  # of x axis
    y_max: float = 1.05,
    figsize: tuple = (12, 6),
    plot_type: str = "box",  # "bar", "box"
    metrics_include: list = ["F1", "Precision", "Recall"],
):
    if not plot_order:
        plot_order = (
            pd.concat(
                [
                    df_results_long.groupby("Feature").Value.median().rename("median"),
                    df_results_long.groupby("Feature").Value.std().rename("std"),
                ],
                axis=1,
            )
            .sort_values(["median", "std"], ascending=[True, False])
            .index
        )

    df_results_long_plt = df_results_long[df_results_long.Metric.isin(metrics_include)]
    df_results_long_plt = df_results_long_plt.assign(
        Metric=df_results_long_plt.Metric.str.replace("_", " ").str.title()
    )

    multiple_metrics = len(metrics_include) > 1
    if not multiple_metrics:
        metric_name = metrics_include[0].replace("_", " ").title()
        df_results_long_plt = df_results_long_plt.drop("Metric", axis=1).rename(
            columns={"Value": metric_name}
        )
    plt.figure(figsize=figsize, dpi=300)
    match plot_type:
        case "box":
            plot_func = sns.boxplot
        case "bar":
            plot_func = sns.barplot
        case x:
            raise ValueError(f"unknown plot type {x}")
    plot_func(
        df_results_long_plt,
        x="Feature",
        y="Value" if multiple_metrics else metric_name,
        hue="Metric" if multiple_metrics else None,
        order=plot_order,
    )
    plt.xticks(rotation=90)
    plt.ylim((0, y_max))
    plt.grid(True, alpha=0.5)
    plt.yticks(np.arange(0, 1.1, 0.1))
    metrics_str = "_".join(metrics_include).replace(" ", "-")
    plt.savefig(
        output_folder_path + test_name + "_" + metrics_str, bbox_inches="tight", dpi=300
    )


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

    outliers = np.array([find_outliers(ml_dataset.X) for ml_dataset in ml_datasets])

    df_outliers = (
        pd.Series(index=ml_datasets[0].sample_names, data=outliers.sum(axis=0))
        .sort_values(ascending=False)
        .to_frame(name="outlier_count")
        .join(dataset_full[0].protein_names)
        .join(dataset_full[1])
    )

    df_outliers = df_outliers[df_outliers.outlier_count >= len(ml_datasets) * threshold]
    # potential outliers to manually check TODO useful?
    return df_outliers
