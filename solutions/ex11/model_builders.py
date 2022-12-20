from typing import Dict, Any, Tuple

from tensorflow import keras
from keras.layers import Dense, Input
from numpy.typing import NDArray


def model_builder_11_2(hparams: Dict[str, Any], input_shape: NDArray,
                       metric: str) -> keras.Sequential:
    """Builds a simple fully connected dnn from the specified hparams to be used in exercise 11.2

    Args:
        hparams (Dict[str, Any]): Set of hyperparameters values to use to build the model.
        input_shape (NDArray): Input shape.
        metric (str): Metric to use to evaluate the model's performance.
    Returns:
        keras.Sequential: Trained model.
    """
    model = keras.Sequential([Input(input_shape)])
    for _ in range(hparams["n_layers"]):
        model.add(
            Dense(
                hparams["n_units"],
                activation=hparams["activation_fn"],
            )
        )
    # Output
    model.add(Dense(1))
    optimizer = keras.optimizers.deserialize(
        {"class_name": hparams["optimizer"],
            "config": {"lr": hparams["lr"]}}
    )
    model.compile(optimizer=optimizer,
                  loss=hparams["loss_fn"], metrics=[metric])
    return model


def model_builder_11_3(hparams: Dict[str, Any], input_shape: NDArray, metric: str) -> keras.Sequential:
    """Builds a simple fully connected dnn from the specified hparams to be used in exercise 11.3

    Args:
        hparams (Dict[str, Any]): Set of hyperparameters values to use to build the model.
        input_shape (NDArray): Input shape.
        metric (str): Metric to use to evaluate the model's performance.
    Returns:
        keras.Sequential: Trained model.
    """
    model = keras.Sequential([Input(input_shape)])
    for _ in range(hparams["n_layers"] - 1):
        model.add(
            Dense(
                hparams["n_units"],
                activation=hparams["activation_fn"],
            )
        )
    model.add(Dense(1))
    optimizer = keras.optimizers.deserialize(
        {"class_name": hparams["optimizer"], "config": {"lr": hparams["lr"]}}
    )
    model.compile(optimizer=optimizer, loss="mae",
                  metrics=[metric])
    return model
