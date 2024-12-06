"""
Parallel Computation Module using Ray

This module provides tools to efficiently parallelize and manage computational tasks
using Ray, a high-performance distributed execution framework. It is designed to
facilitate the distribution of workloads across multiple CPU cores or machines,
thereby accelerating processing times for large-scale data or compute-intensive
operations.

## Usage

1. **Define a Task Class**

    Subclass `RayTaskGeneral` and implement the `process_item` method with the desired computation.

    ```python
    class MyTask(RayTaskGeneral):
        def process_item(self, item):
            # Implement the computation
            return item * 2  # Example operation
    ```

2. **Initialize Ray and RayManager**

    ```python
    import ray
    from your_module import RayManager, MyTask

    ray.init()  # Initialize Ray

    manager = RayManager(task_class=MyTask, n_cores=4)  # Use 4 cores
    ```

3. **Submit Tasks**

    ```python
    items_to_process = [1, 2, 3, 4, 5]

    def save_results(results, **kwargs):
        # Implement saving logic, e.g., write to a file or database
        print("Saving results:", results)

    manager.submit_tasks(
        items=items_to_process,
        chunk_size=2,
        save_func=save_results,
        save_kwargs={'destination': 'output.txt'},
        save_interval=2
    )
    ```

4. **Retrieve Results**

    ```python
    all_results = manager.get_results()
    print("All Results:", all_results)
    ```

## Notes

The `RayManager` handles task submission and result collection efficiently by
maintaining a pool of worker futures up to the specified number of cores.
The `save_func` allows for intermediate saving of results, which is useful for
long-running tasks to prevent data loss and manage memory usage.
Exception handling is incorporated to ensure that individual task failures do not
halt the entire processing pipeline. Errors are logged and returned as part of
the results for further inspection.


"""

from typing import Any, Generator
from abc import ABC, abstractmethod

from collections.abc import Callable

import ray
from loguru import logger


@ray.remote
def ray_worker(task_class: Callable[[], Any], item: Any) -> Any:
    """
    Remote function to process a single item using the provided task class.

    Args:
        task_class (Callable[[], Any]): A callable that returns an instance with a `run` method.
        item (Any): The input data required for the calculation.

    Returns:
        Any: The result of the `run` method or an error tuple.
    """
    task_instance = task_class()
    result = task_instance.run(item)
    return result


class RayManager:
    """
    A manager class for handling task submissions and result collection using Ray's Task Parallelism.
    """

    def __init__(self, task_class: Callable[[], Any], n_cores: int = -1) -> None:
        """Initializes the RayManager.

        Args:
            task_class: A callable that returns an instance with a `run` method for
                processing each item.
            n_cores: Number of parallel tasks to run. If <= 0, uses all available CPUs.
        """
        self.task_class = task_class
        """
        A callable that returns an instance with a `run` method for
        processing each item.
        """

        self.n_cores = (
            n_cores if n_cores > 0 else int(ray.available_resources().get("CPU", 1))
        )
        """
        The number of parallel tasks to run. If set to `-1` or any value less than or
        equal to `0`, all available CPU cores are utilized.
        """

        self.futures: list[ray.ObjectRef] = []
        """
        A list of Ray object references representing the currently submitted but
        not yet completed tasks. This manages the pool of active workers.
        """

        self.results: list[Any] = []
        """
        A list that stores the results of all completed tasks. It aggregates the
        output returned by each worker.
        """

        self.save_func: Callable[[Any, Any], None] | None = None
        """
        An optional callable function that takes a batch of results and performs a
        save operation. This can be used to persist intermediate results to
        disk, a database, or any other storage medium. If set to `None`, results are
        not saved automatically.
        """

        self.save_interval: int = 1
        """
        The number of results to accumulate before invoking `save_func`.
        When the number of collected results reaches this interval,
        `save_func` is called to handle the batch of results.
        """

        self.save_kwargs: dict[str, Any] = dict()
        """
        A dictionary of additional keyword arguments to pass to
        `save_func` when it is called. This allows for flexible configuration of the
        save operation, such as specifying file paths,
        database connections, or other parameters required by `save_func`.
        """

    def task_generator(
        self, items: list[Any], chunk_size: int
    ) -> Generator[Any, None, None]:
        """Generator that yields individual items from chunks.

        Args:
            items: A list of items to process.
            chunk_size: Number of items per chunk.

        Yields:
            Individual items to be processed.
        """
        for start in range(0, len(items), chunk_size):
            end = min(start + chunk_size, len(items))
            chunk = items[start:end]
            logger.info(f"Yielding tasks {start} to {end - 1}")
            for item in chunk:
                yield item

    def submit_tasks(
        self,
        items: list[Any],
        chunk_size: int = 100,
        save_func: Callable[[Any, Any], None] | None = None,
        save_kwargs: dict[str, Any] = dict(),
        save_interval: int = 100,
    ) -> None:
        """Submits tasks using a generator and manages workers up to n_cores.

        Args:
            items: A list of items to process.
            chunk_size: Number of items per chunk.
            save_func: A callable that takes a list of results and saves them.
            save_interval: The number of results after which to invoke save_func.
        """
        self.save_func = save_func
        self.save_kwargs = save_kwargs
        self.save_interval = save_interval

        task_gen = self.task_generator(items, chunk_size)

        results = []
        for item in task_gen:
            if len(self.futures) >= self.n_cores:
                # Wait for any worker to finish
                done_futures, self.futures = ray.wait(self.futures, num_returns=1)
                result = ray.get(done_futures[0])
                results.append(result)

                # Save results if save_func is set and interval is reached
                if self.save_func and len(results) >= self.save_interval:
                    self.save_func(results[: self.save_interval], **self.save_kwargs)
                    self.results.extend(results)
                    results = results[self.save_interval :]

            # Submit new task
            future = ray_worker.remote(self.task_class, item)
            self.futures.append(future)

        # Collect remaining futures
        while self.futures:
            done_futures, self.futures = ray.wait(self.futures, num_returns=1)
            result = ray.get(done_futures[0])
            results.append(result)

            # Save results if save_func is set and interval is reached
            if self.save_func and len(results) >= self.save_interval:
                self.save_func(results[: self.save_interval], **self.save_kwargs)
                self.results.extend(results)
                results = results[self.save_interval :]

        # Save any remaining results
        if self.save_func and results:
            self.save_func(results, **self.save_kwargs)
        self.results.extend(results)

    def get_results(self) -> list[Any]:
        """
        Retrieves all collected results.

        Returns:
            A list of results from all completed tasks.
        """
        return self.results


class RayTaskGeneral(ABC):
    """A parent class of others that governs what calculations are run on each
    task."""

    def __init__(self):
        self.results: list[Any] = []

    def run(self, item: Any) -> Any:
        """Processes a single item.

        This method wraps the `process_item` method with a try-except block.

        Args:
            item: The input data required for the calculation.

        Returns:
            The result of the calculation or an error tuple.
        """
        try:
            result = self.process_item(item)
            return result
        except Exception as e:
            logger.exception(f"Error processing item {item}: {e}")
            return ("error", str(e))

    @abstractmethod
    def process_item(self, item: Any) -> Any:
        """The definition that computes or processes a single item.

        This method should be implemented by subclasses to define the specific
        processing logic for each item.

        Args:
            item: The input data required for the calculation.

        Returns:
            Any: The result of processing the item.
        """
