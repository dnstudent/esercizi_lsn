from tensorflow import keras
import tensorflow as tf
from tensorboard.plugins.hparams import api as hp
from numpy.typing import NDArray
import pandas as pd
import numpy as np
from typing import Tuple, Dict, Any, Callable, List, Union, Optional
DataSet = Tuple[NDArray, NDArray]


def train_model(model_builder: Callable,
                hparams: Dict[str, Any],
                metric: str,
                train_data: DataSet,
                valid_data: Optional[DataSet],
                n_epochs: int,
                batch_size: int,
                early_stopping: Optional[str],
                callbacks: List[keras.callbacks.Callback] = [],
                **es_kwargs):
    if type(hparams) is pd.Series:
        hparams = hparams.to_dict()
    x_train, y_train = train_data
    # building the model
    model = model_builder(hparams, np.shape(x_train)[1:], metric)
    validation_kws = {"validation_data": None}
    # defining the callbacks:
    # EarlyStopping
    if early_stopping is not None:
        # Defining which metric will be monitored by the callback
        if early_stopping == "train":
            monitor = "loss"
        elif early_stopping == "valid":
            assert (valid_data is not None)
            monitor = f"val_loss"
            # Defining additional parameters to the fit function: early stopping needs a validation dataset
            validation_kws["validation_data"] = valid_data
        else:
            raise NotImplementedError(
                f"Argument {early_stopping} is not supported for early_stopping. Only ones are 'train' and 'valid'")
        callbacks.append(
            keras.callbacks.EarlyStopping(
                monitor, mode="min", **es_kwargs)
        )

    # fitting the model
    history = model.fit(
        x_train,
        y_train,
        epochs=n_epochs,
        batch_size=batch_size,
        verbose=0,
        shuffle=True,
        callbacks=callbacks,
        **validation_kws
    )
    return model, history


def model_initializer(model_builder: Callable[[Dict[str, Any], NDArray, str], keras.Model],
                      hparams: Union[Dict[str, Any], pd.Series],
                      metric: str,
                      train_data: DataSet,
                      valid_data: DataSet,
                      n_epochs: int,
                      batch_size: int,
                      early_stopping: Optional[str],
                      callbacks: List[keras.callbacks.Callback] = [],
                      **es_kwargs):
    model, history = train_model(model_builder, hparams, metric,
                                 train_data, valid_data, n_epochs, batch_size, early_stopping, callbacks, **es_kwargs)
    # The model is evaluated only if it hasn't been evaluated during the fit
    if early_stopping != "valid":
        valid_loss, valid_metric = model.evaluate(
            valid_data[0], valid_data[1], verbose=0)
    else:
        if "restore_best_weights" in es_kwargs and es_kwargs["restore_best_weights"]:
            idx = np.argmin(history.history[f"val_{metric}"])
        else:
            idx = -1
        valid_loss = history.history["val_loss"][idx]
        valid_metric = history.history[f"val_{metric}"][idx]
    return model, {"history": history.history, "val_loss": valid_loss, f"val_{metric}": valid_metric}


def eval_params(model_builder: Callable,
                hparams: Dict[str, Any],
                metric: str,
                train_data: DataSet,
                valid_data: DataSet,
                n_epochs: int,
                batch_size: int,
                early_stopping: Optional[str],
                callbacks: List[keras.callbacks.Callback] = [], **es_kwargs) -> Tuple[keras.Model, float, float, float, float]:
    """Evaluates the model built using the given hyperparameters using the given metric.

    Args:
        model_builder (Callable): A model builder taking at least hparams and a metric in input and returning a Keras model and its validation history.
        hparams (Dict[str, Any]): A set of hyperparameters.
        metric (string): A metric to evaluate the model.
        **kwargs: Additional arguments passed to model_initializer.

    Returns:
        Tuple[float, float]: Validation loss and metric
    """
    model, score = model_initializer(model_builder, hparams, metric, train_data,
                                     valid_data, n_epochs, batch_size, early_stopping, callbacks, **es_kwargs)
    return model, {"loss": score["history"]["loss"][-1], "val_loss": score["val_loss"], metric: score["history"][metric][-1], f"val_{metric}": score[f"val_{metric}"]}


def evaluate_and_log(model_builder: Callable,
                     hparams: Dict[str, Any],
                     metric: str,
                     logdir: str,
                     train_data: DataSet,
                     valid_data: DataSet,
                     n_epochs: int,
                     early_stopping: Optional[str],
                     callbacks: List[keras.callbacks.Callback] = [], **es_kwargs):
    """Evaluates the given set of hyperparameters and logs the results in Tensorboard

    Args:
        model_builder (Callable): A model builder taking at least hparams and a metric in input and returning a Keras model and its validation history.
        hparams (Dict[str, Any]): A set of hyperparameters.
        metric (string): A metric to evaluate the model.
        logdir (str): The path in which logs will be stored.
        **kwargs: Additional arguments passed to model_initializer.

    Returns:
        Tuple[float, float]: Validation loss and metric
    """
    with tf.summary.create_file_writer(logdir).as_default():
        hp.hparams(hparams)  # record the values used in this trial
        _, scores = eval_params(
            model_builder, hparams, metric, train_data, valid_data, n_epochs, early_stopping, callbacks, es_kwargs)
        tf.summary.scalar(metric, scores[metric], step=1)
    return scores
