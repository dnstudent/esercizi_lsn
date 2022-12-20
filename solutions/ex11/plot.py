from typing import Optional, Callable, Tuple, List, Dict

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from tensorflow import keras

Dataset = Tuple[NDArray, NDArray]

# 11.2 section


def show_metrics(metric, search_results_df, axes, color_by):
    for dataset, ax in zip(["training_", "inside_", "outside_"], axes):
        sns.swarmplot(search_results_df,
                      x="model_type",
                      y=f"{dataset}{metric}",
                      hue=color_by,
                      legend=(ax == axes[0]),
                      ax=ax)
        ax.grid(which="both")


def make_base_plotter(xs: NDArray, training_dataset: Dataset, test_datasets: Dict[str, Dataset], target: Callable):
    def show_base(ax: plt.Axes):
        """Plots the target function evaluated in a given domain, the training and test datasets.

        Args:
            ax (plt.Axes): Matplotlib Axes where everything will be plotted.
        """
        handles = []
        handles.append(ax.plot(xs, target(xs), c="orange", label="f(x)")[0])
        handles.append(ax.scatter(
            training_dataset[0], training_dataset[1], s=5, c="r", label="training dataset"))
        handles.append(ax.scatter(
            test_datasets["inside"][0], test_datasets["inside"][1], c="g", s=0.5, label="validation datasets"))
        ax.scatter(test_datasets["outside"][0],
                   test_datasets["outside"][1], c="g", s=0.5)
        # Returns the handles used to show a legend
        return handles
    return show_base


def make_model_plotter(xs: NDArray, models: List[keras.Model]):
    def show_model_predictions(ax: plt.Axes, model_number: int, color: str = "b", label: Optional[str] = None):
        """Plots a model's predictions.

        Args:
            ax (plt.Axes): The Axes where the plot will be drawn.
            model_number (int): Model index in the "models" list.
            color (str, optional): Matplotlib color for the plot. Defaults to "b".
            label (Optional[str], optional): Legend label for the plot. Defaults to None.
        """
        m = models[model_number]
        if label is None:
            ax.plot(xs, m(xs), c=color)
            return []
        return ax.plot(xs, m(xs), c=color, label=label)
    return show_model_predictions


def make_panels_plotter(xs: NDArray, models: keras.Model, training_dataset: Dataset, test_datasets: Dict[str, Dataset], target: Callable[[NDArray], NDArray]) -> Callable:
    """Generates the panel plotting function.

    Args:
        xs (NDArray): Coordinates where the target function will be evaluated.
        models (keras.Model): List of models.
        training_dataset (Dataset): Training dataset.
        test_datasets (Dict[str, Dataset]): Test datasets ("inside" and "outside").
        target (Callable[[NDArray], NDArray]): Target function.

    Returns:
        Callable: Panel plotting function.
    """
    show_base = make_base_plotter(xs, training_dataset, test_datasets, target)
    show_model_predictions = make_model_plotter(xs, models)

    def plot_panels(model_numbers: Tuple[List[int], List[int]], title: str, colors: List[str] = ["b", "c"]) -> List[plt.Axes]:
        """Draws a panel consisting of two plots: on the left the one with "simple" models results, on the right with "complex" models results.
        For each model kind the specified models are taken from the list of models and their predictions on the set "xs" are drawn.

        Args:
            model_numbers (Tuple[List[int], List[int]]): Model indices to select models from the "models" list.
            title (str): Title to give to the panel.
            colors (List[str], optional): List of styles used to plot models predictions. Must be one for each selected model. Defaults to ["b", "c"].

        Returns:
            List[plt.Axes]: The panel's Axes
        """

        _, axes = plt.subplots(1, 2, figsize=(
            10, 5), constrained_layout=True, sharey=True)
        # Iterating over the panel's axes
        for model_kind, panel_numbers, ax in zip(["Simple", "Complex"], model_numbers, axes):
            # Plotting everything and collecting legend handles
            handles = (show_base(ax) +
                       sum([show_model_predictions(ax, number, color, label=f"model_{number}")
                            for number, color in zip(panel_numbers, colors)], start=[]))
            ax.legend(handles=handles)
            _ = ax.set_title(f"{model_kind} models {title}")
        return axes
    return plot_panels


# 11.3 section
def plot_target(xs: NDArray, ys: NDArray, target: Callable, n_levels: int, ax: plt.Axes):
    """Plots a contour plot of the target function on an ax

    Args:
        xs, ys (NDArray): List of x, y coordinates where the target function will be computed.
        target (Callable): Target function
        ax (plt.Axes): Matplotlib axes
    """
    zs_target = target(np.stack([xs, ys], axis=-1))
    c1 = ax.contourf(
        xs,
        ys,
        zs_target,
        levels=np.linspace(-1.5, 1.5, n_levels, endpoint=True),
        cmap="bwr",
    )
    ax.set_title(r"f(x,y)")
    ax.figure.colorbar(c1, ax=ax)
    return zs_target


def plot_model_predictions(xs: NDArray, ys: NDArray, model: keras.Model, n_levels: int, ax: plt.Axes) -> NDArray:
    """Computes and plots a model's predictions using a countour plot.

    Args:
        xs (NDArray): List of x coordinates where the target function will be computed.
        ys (NDArray): List of y coordinates where the target function will be computed.
        model (keras.Model): Model.
        n_levels (int): Number of levels in the countour plot.
        ax (plt.Axes): Matplotlib Axes.

    Returns:
        NDArray: Model's predictions.
    """
    zs_interp = model.predict(
        np.stack([xs, ys], axis=-1).reshape((-1, 2)), verbose=0, batch_size=None
    ).reshape(xs.shape)
    c2 = ax.contourf(
        xs,
        ys,
        np.clip(zs_interp, -1.5, 1.5),
        levels=np.linspace(-1.5, 1.5, n_levels, endpoint=True),
        cmap="bwr",
    )
    ax.figure.colorbar(c2, ax=ax)
    return zs_interp


def plot_difference(xs: NDArray,
                    ys: NDArray,
                    zs_target: NDArray,
                    interp: NDArray,
                    n_levels: int,
                    train_data: Tuple[NDArray, NDArray],
                    valid_data: Tuple[NDArray, NDArray],
                    ax: plt.Axes):
    """Plots the difference between target and predicted values using a contour plot.

    Args:
        xs (NDArray): List of x coordinates where the target function will be computed.
        ys (NDArray): List of y coordinates where the target function will be computed.
        zs_target (NDArray): List of z coordinates evaluated by the target function.
        interp (NDArray): List of z coordinates estimated by a model.
        n_levels (int): Number of levels in the contour plot.
        train_data (Tuple[NDArray, NDArray]): Training dataset.
        valid_data (Tuple[NDArray, NDArray]): Validation dataset.
        ax (plt.Axes): Axes where the plot will be drawn.
    """
    diff = zs_target - interp
    clip = np.max(np.abs(diff))
    c3 = ax.contourf(
        xs,
        ys,
        diff,
        levels=np.linspace(-clip, clip, n_levels, endpoint=True),
        cmap="bwr",
    )
    p1 = ax.scatter(
        train_data[0][:, 0], train_data[0][:, 1], s=0.2, c="g", label="training_data")
    p2 = ax.scatter(valid_data[0][:, 0], valid_data[0]
                    [:, 1], s=0.2, c="black", label="validation_data")
    ax.set_title("Difference + training points")
    ax.legend(handles=[p1, p2], loc="upper left")
    _ = ax.figure.colorbar(c3, ax=ax)


def plot_panel(xs, ys, zs_target, model, model_type, n_levels, train_data, valid_data, axes: List[plt.Axes]):
    """Plots both the model predictions and the difference with the target on two separate axes. Check the documentations from plot_model_predictions and plt_difference
    """
    zs_interp = plot_model_predictions(xs, ys, model, n_levels, axes[0])
    axes[0].set_title(model_type)
    plot_difference(xs, ys, zs_target, zs_interp, n_levels,
                    train_data, valid_data, axes[1])
