from sklearn.preprocessing import StandardScaler
from keras.layers import Dense
from keras.models import Sequential
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import KFold, cross_val_score
from sklearn.pipeline import Pipeline
import hyperopt
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials


def bl(space):
    model = Sequential()
    model.add(Dense(3, input_dim=3, kernel_initializer="normal",
                    activation="relu"))
    model.add(Dense(3, kernel_initializer="normal"))
    if space["n_hlayers"] > 1:
        if space["n_hlayers"] == 2:
            model.add(Dense(3, kernel_initializer="normal"))
        elif space["n_hlayers"] == 3:
            model.add(Dense(3, kernel_initializer="normal"))
            model.add(Dense(3, kernel_initializer="normal"))
        elif space["n_hlayers"] == 4:
            model.add(Dense(3, kernel_initializer="normal"))
            model.add(Dense(3, kernel_initializer="normal"))
            model.add(Dense(3, kernel_initializer="normal"))
        else:
            raise Exception("Should be 1, 2, 3 or 4 HL")
    model.add(Dense(1, kernel_initializer="normal"))
    model.compile(loss="mean_squared_error", optimizer="adam")

    return model


estimators = []
estimators.append(("stand", StandardScaler()))
estimators.append(("mlp", KerasRegressor(build_fn=bl, epochs=100, batch_size=10)))
pipeline = Pipeline(estimators)
kfold = KFold(n_splits=5, random_state=42)