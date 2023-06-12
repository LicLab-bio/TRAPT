import random
import sys

import numpy as np
import pandas as pd
import tensorflow as tf
from keras import Input, Model
from keras import backend as K
from keras.callbacks import EarlyStopping
from keras.layers import Dense, Dropout
from keras.models import Model
from keras.optimizers import Adam
from keras.regularizers import Regularizer

def seed_tensorflow(seed=2023):
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)


def custom_sigmoid(x, t=1):
    return 1 / (1 + K.exp(-x / t))


class SparseGroupLasso(Regularizer):
    """
    https://keras-ssg-lasso.readthedocs.io/en/latest/
    argmin(w) (1/2n * ||y - Xw||^2 + λ1 * ||w||1 + λ2 * ∑g∈G √pg * ||w_g||2)
    """

    def __init__(self, l1=0.01, l2=0.01, groups=None):
        self.l1 = l1
        self.l2 = l2
        groups = np.array(groups).astype("int32")
        self.p_l = K.variable(np.sqrt(np.bincount(groups)).reshape((1, -1)))
        self.groups = K.variable(groups, "int32")

    def __call__(self, x):
        regularization = 0.0
        if self.l1:
            regularization += self.l1 * K.sum(K.abs(x))
        if self.l2:
            regularization += self.l2 * K.sum(
                K.dot(
                    self.p_l,
                    K.sqrt(
                        tf.math.unsorted_segment_sum(
                            K.square(x), self.groups, num_segments=self.p_l.shape[1]
                        )
                    ),
                )
            )
        return regularization

    def get_config(self):
        return {
            "name": self.__class__.__name__,
            "l1": self.l1,
            "l2": self.l2,
            "groups": self.groups,
        }


class FeatureSelection:
    def __init__(self, args, data_ad, type):
        self.args = args
        self.data_ad = data_ad
        self.type = type
        self.epochs = 16
        self.batch_size = 32
        self.learning_rate = 0.001
        self.group_size = 10
        self.background_genes = 6000
        self.top = 10
        self.shuffle = True
        self.is_val = False
        self.early_stopping = EarlyStopping(patience=4)

    def TSFS(self, X, T):
        seed_tensorflow()
        pos_weight = float(len(T) - T.sum()) / T.sum()
        sample_weight = np.ones(len(T))
        sample_weight[T.astype(bool)] = pos_weight
        input = Input(X.shape[1:])
        y = Dense(2, activation="relu")(input)
        feature_extraction = Model(input, y)
        t = feature_extraction(input)
        output = Dense(1)(t)
        output = custom_sigmoid(output, 1)
        teacher = Model(input, output)
        teacher.compile(
            optimizer=Adam(self.learning_rate), loss="mse", weighted_metrics=[]
        )
        if self.is_val:
            teacher.fit(
                X,
                T,
                epochs=self.epochs,
                batch_size=self.batch_size,
                callbacks=[self.early_stopping],
                validation_split=0.15,
                sample_weight=sample_weight,
                shuffle=self.shuffle,
                verbose=1,
            )
        else:
            teacher.fit(
                X,
                T,
                epochs=self.epochs,
                batch_size=self.batch_size,
                sample_weight=sample_weight,
                shuffle=self.shuffle,
                verbose=0,
            )
        Y = feature_extraction.predict(X)
        seed_tensorflow()
        input = Input((X.shape[-1],))
        groups = self.groups
        output = Dense(
            Y.shape[-1] * 10,
            activation="relu",
            kernel_regularizer=SparseGroupLasso(groups=groups),
        )(input)
        o_1 = Dense(Y.shape[-1], activation="relu")(output)
        student = Model(input, o_1)
        student.compile(optimizer=Adam(self.learning_rate), loss="mse")
        student.fit(
            X,
            Y,
            epochs=self.epochs,
            batch_size=self.batch_size,
            sample_weight=sample_weight,
            shuffle=self.shuffle,
            verbose=0,
        )
        weights = student.layers[1].get_weights()[0]
        weights = np.nan_to_num(weights, nan=0)
        weights = np.sum(np.square(weights), 1)
        return np.argsort(-weights)

    def train(self, X, y):
        seed_tensorflow()
        pos_weight = float(len(y) - y.sum()) / y.sum()
        sample_weight = np.ones(len(y))
        sample_weight[y.astype(bool)] = pos_weight
        input = Input(X.shape[1:])
        input = Dropout(0.1)(input)
        output = Dense(2, activation="relu")(input)
        output = Dense(1)(output)
        output = custom_sigmoid(output, 1)
        model = Model(input, output)
        model.compile(
            optimizer=Adam(self.learning_rate), loss="mse", weighted_metrics=[]
        )
        if self.is_val:
            model.fit(
                X,
                y,
                epochs=self.epochs,
                callbacks=[self.early_stopping],
                validation_split=0.15,
                batch_size=self.batch_size,
                sample_weight=sample_weight,
                shuffle=self.shuffle,
                verbose=1,
            )
        else:
            model.fit(
                X,
                y,
                epochs=self.epochs,
                batch_size=self.batch_size,
                sample_weight=sample_weight,
                shuffle=self.shuffle,
                verbose=0,
            )
        return model

    def get_corr(self, v1: list, v2: list):
        return np.array([np.corrcoef(v1, v2[i])[0, 1] for i in range(len(v2))])

    def sort_by_group(self, vec):
        vec = self.get_corr(vec, self.data_ad.X)
        self.index = np.argsort((vec.min() - vec))
        self.data_ad = self.data_ad[self.index]
        self.groups = [i // self.group_size for i in range(len(self.index))]

    def run(self):
        geneset = pd.read_csv(self.args.input, header=None)[0]
        genes = self.data_ad.var_names.values
        y = np.in1d(genes, geneset).astype(int)
        self.sort_by_group(y)
        X = self.data_ad.X.T
        true_index = np.where(y.astype(bool))[0]
        np.random.seed(2023)
        false_index = np.random.choice(
            np.where(~y.astype(bool))[0], self.background_genes, replace=False
        )
        index = np.append(true_index, false_index)
        feature_indexs = self.TSFS(X[index], y[index])
        feature_indexs = feature_indexs[: self.top]
        features = X[:, feature_indexs]
        model = self.train(features[index], y[index])
        gene_vec = model.predict(features)
        gene_vec = pd.DataFrame(gene_vec, index=self.data_ad.var_names)
        return gene_vec
