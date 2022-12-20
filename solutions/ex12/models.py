import re
import os
from pathlib import Path

from tensorflow import keras
from keras.layers import Conv2D, Dense, Dropout, Flatten, Input, MaxPooling2D


def make_dnn(name, input_shape=(28 * 28,), n_classes=10, **kwargs) -> keras.Sequential:
    """Builds a simple feedforward NN with two fully connected layers and a dropout layer.
    The model computes logits, not a probability distribution.

    Returns:
        keras.Sequential: The model itself
    """
    return keras.Sequential(
        [
            Input(input_shape),
            Dense(512, activation="relu"),
            Dense(128, activation="relu"),
            Dropout(0.1),
            Dense(n_classes),
        ],
        name=name,
        **kwargs,
    )


def make_small_cnn(
    name, input_shape=(28, 28, 1), n_classes=10, **kwargs
) -> keras.Sequential:
    """Builds a simple convolutional NN with dropout layer.
    The model computes logits, not a probability distribution.
    It is aimed at processing mnist digits.

    Returns:
        keras.Sequential: The model itself
    """
    return keras.Sequential(
        [
            Input(input_shape),
            Conv2D(
                32,
                kernel_size=3,
                activation="relu",
                padding="valid",
                data_format="channels_last",
            ),
            MaxPooling2D(2),
            Flatten(),
            Dense(128, activation="relu"),
            Dropout(0.1),
            Dense(n_classes),
        ],
        name=name,
        **kwargs,
    )


def make_big_cnn(
    name, input_shape=(28, 28, 1), n_classes=10, **kwargs
) -> keras.Sequential:
    """Builds a convolutional NN with dropout layer.
    The model computes logits, not a probability distribution.
    It is aimed at processing mnist digits.

    Returns:
        keras.Sequential: The model itself
    """
    return keras.Sequential(
        [
            Input(input_shape),
            Conv2D(
                32,
                kernel_size=5,
                activation="relu",
                padding="valid",
                data_format="channels_last",
            ),
            MaxPooling2D(2),
            Conv2D(
                64,
                kernel_size=5,
                activation="relu",
                padding="valid",
            ),
            MaxPooling2D(2),
            Dropout(0.4),
            Flatten(),
            Dense(128, activation="relu"),
            Dropout(0.4),
            Dense(n_classes),
        ],
        name=name,
        **kwargs,
    )


def storage_path(model: keras.Model) -> Path:
    """Returns the path which should store a model's weights in my arrangement

    Args:
        model (keras.Model): A Keras model.
    """
    os.makedirs("models", exist_ok=True)
    return Path(os.path.join("models", f"{model.name}"))


def _extract_acc(savepath):
    data = re.search("\d+-([\d\.]+).hdf5$", str(savepath))
    return float(data.groups()[0])


def checkpoint_exists(model: keras.Model) -> bool:
    """Checks whether a checkpoint of the model's weights exists

    Args:
        model (keras.Model): A Keras Model.
    """
    models_path = storage_path(model)
    return len(list(models_path.glob("**/*.hdf5"))) > 0


def load_best_weights(model: keras.Model, verbose=0) -> bool:
    """Checks if weights for a model by the name of `name` exist. If so they get loaded.

    Args:
        model (keras.Model): The model on which the weights will be loaded.
        verbose (int, optional): Level of output verbosity. Defaults to 0.

    Returns:
        bool: Whether a checkpoint was loaded.
    """
    models_path = storage_path(model)
    files = list(models_path.glob("**/*.hdf5"))
    if len(files) > 0:
        ckpt = max(files, key=_extract_acc)
        try:
            model.load_weights(ckpt)
            # data = re.search("(\d+)-([\d\.]+).hdf5$", str(ckpt))
            # epoch = int(data.groups()[0])
            # val_accuracy = float(data.groups()[1])
            found = True
        except:
            found = False
    else:
        found = False
    if verbose > 0:
        if found:
            print(f"{model.name} was found")
        else:
            print(f"{model.name} was not found")
    return found
