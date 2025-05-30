"""
@author: adenger
"""

from dataclasses import dataclass
import numpy as np
from sklearn.preprocessing import LabelEncoder

@dataclass
class MLDataset:
    name: str
    X: np.array
    y: np.array
    feature_names: np.array
    sample_names: np.array
    label_encoder: LabelEncoder

    def __str__(self):
        n_classes = self.label_encoder.transform(self.label_encoder.classes_)
        n_labels = self.label_encoder.classes_
        n_features = len(self.feature_names)
        n_samples = len(self.sample_names)
        return f"Name: {self.name}, Features: {n_features}, Samples: {n_samples}, Classes: {n_classes}, Labels: {n_labels})"

    def __repr__(self):
        return self.__str__()