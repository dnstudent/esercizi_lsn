from typing import Union, Tuple, List, Dict, Optional
import os
from pathlib import Path

from sklearn import metrics
import matplotlib.pyplot as plt
from tensorflow import keras
import numpy as np

from global_utils import results_dir

RESULTS_PATH = results_dir("12")


def history_path(model: keras.Model) -> Path:
    os.makedirs(RESULTS_PATH / "metrics", exist_ok=True)
    return RESULTS_PATH / "metrics" / f"{model.name}.npz"


def save_history(model: keras.Model, history: keras.callbacks.History, overwrite: bool):
    """Stores a model's training history

    Args:
        model (keras.Model): A Keras Model.
        history (keras.callbacks.History): The training history.
        overwrite (bool): Should previous history be overwritten?
    """
    path = history_path(model)
    if overwrite or not os.path.exists(path):
        np.savez(path, **history.history)


def load_history(model: keras.Model) -> Optional[keras.callbacks.History]:
    """Load a model's training history

    Args:
        model (keras.Model): A Keras Model.

    Returns:
        keras.callbacks.History: The model's training history.
    """
    path = history_path(model)
    try:
        history = keras.callbacks.History()
        with np.load(path) as f:
            history.history = dict(f)
        history.epoch = list(range(len(history.history["loss"])))
    except:
        history = None
    return history


def optimize_model(
    model: keras.Model,
    data: Tuple[Tuple[np.array, np.array], Tuple[np.array, np.array]],
    optimizer: Union[str, keras.optimizers.Optimizer],
    batch_size,
    starting_val_accuracy=0.0,
    verbose=0,
    max_epochs=15,
) -> keras.callbacks.History:
    """Evaluates the performance of a model fitted with the specified optimizer.
    It implements early stopping.

    Args:
        data: A tuple containing ((x_train, y_train), (x_test, y_test)) as numpy arrays.
        model (keras.Model): A Keras model. It is required to have a name.
        optimizer (Union[str, keras.optimizers.Optimizer]): A valid Keras optimizer.
        batch_size (int): The batch size used for training.
        verbose (int, optional): Verbosity of model.fit. Defaults to 0.
        max_epochs (int, optional): Maximum number of epochs. Defaults to 15.
    Returns:
        keras.callbacks.History: Fit history, complete with per-epoch validation metrics
    """
    (x_train, y_train), (x_test, y_test) = data
    model.compile(
        loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
        optimizer=optimizer,
        metrics=["accuracy"],
    )
    # The model is fit with EarlyStopping on the validation accuracy, preventing it to waste time overfitting
    history = model.fit(
        x_train,
        y_train,
        epochs=max_epochs,
        batch_size=batch_size,
        shuffle=True,
        validation_data=(x_test, y_test),
        callbacks=[
            keras.callbacks.EarlyStopping(
                monitor="val_accuracy",
                # Improvements <= 1e-4 in validation accuracy are considered noise
                min_delta=1e-4,
                mode="max",
                # We wait 3 epochs before declaring overfitting
                patience=10,
                baseline=starting_val_accuracy,
                restore_best_weights=True,
            ),
            #  Checkpointing the model
            keras.callbacks.ModelCheckpoint(
                f'{RESULTS_PATH / "models" / model.name / "{epoch:03d}-{val_accuracy:.4f}.hdf5"}',
                monitor="val_accuracy",
                save_best_only=True,
                mode="max",
                verbose=min(verbose, 1),
            ),
        ],
        verbose=verbose,
    )
    return history


def plot_wrong_predictions(
    fig: plt.Figure,
    model: keras.Model,
    x_test: np.array,
    y_true: np.array,
    rows: int,
    cols: int,
    verbose=0,
):
    y_predicted = np.argmax(model.predict(x_test, verbose=verbose), axis=-1)
    wrong = y_predicted != y_true
    for i, (digit_image, prediction, correct, _) in enumerate(
        zip(x_test[wrong], y_predicted[wrong],
            y_true[wrong], range(rows * cols))
    ):
        ax = fig.add_subplot(rows, cols, i + 1)
        ax.imshow(digit_image.reshape((28, 28, 1)), cmap="gray")
        ax.set_title(f"guess: {prediction} true: {correct}")
        ax.set_axis_off()


def plot_fits(
    training_histories: Dict[str, keras.callbacks.History],
    axes: List[List[plt.Axes]],
):
    """Plots training histories on provided axes.

    Args:
        training_results (Dict[str, keras.callbacks.History]): Training histories.
        axes (List[List[plt.Axes]]): Matplotlib Axes. Must be a 2x2 matrix.
    """
    axes[0, 0].set_title("Training")
    axes[0, 1].set_title("Validation")
    axes[1, 0].set_xlabel("epoch")
    axes[1, 1].set_xlabel("epoch")
    axes[0, 0].set_ylabel("Categorical crossentropy")
    axes[1, 0].set_ylabel("Accuracy")
    for axs, metric in zip(axes, ["loss", "accuracy"]):
        for ax, prefix in zip(axs, ["", "val_"]):
            ax.grid(which="both")
            for name, history in training_histories.items():
                ax.plot(
                    history.epoch,
                    history.history[prefix + metric],
                    label=name,
                )
    axes[1, 1].legend()


def plot_confusion_matrices(
    fig: plt.Figure,
    models: List[keras.Model],
    x_test: np.array,
    y_true: np.array,
    rows: int,
    cols: int,
):
    """Plots confusion matrices for the classifiers in model using provided data

    Args:
        fig (plt.Figure): A Figure object.
        models (List[keras.Model]): Keras Classifiers to be evaluated.
        x_test (np.array): Data to be classified.
        y_test (np.array): True data labels.
        rows (int): Number of rows in the resulting grid.
        cols (int): Number of columns in the resulting grid.
    """
    fig.suptitle("Confusion matrices")
    for i, model in enumerate(models):
        y_predicted = np.argmax(model.predict(x_test, verbose=0), axis=-1)
        ax = fig.add_subplot(rows, cols, i + 1)
        ax.set_title(model.name)
        metrics.ConfusionMatrixDisplay.from_predictions(
            y_true, y_predicted, ax=ax, colorbar=False
        )
