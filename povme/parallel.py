from typing import Any, Generator

from collections.abc import Callable

import ray
from loguru import logger


@ray.remote
def ray_worker(task_class: Callable[[], Any], item: Any) -> Any:
    """
    Remote function to process a single item using the provided task class.

    Args:
        task_class: A callable that returns an instance with a `run` method.
        item: The input data required for the calculation.

    Returns:
        The result of the `run` method or an error tuple.
    """
    try:
        task_instance = task_class()
        result = task_instance.run(item)
        return result
    except Exception as e:
        return ("error", str(e))


class RayManager:
    """
    A manager class for handling task submissions and result collection using Ray's Task Parallelism.
    """

    def __init__(self, task_class: Callable[[], Any], n_cores: int = -1) -> None:
        """Initializes the RayManager.

        Args:
            task_class: A callable that returns an instance with a `run` method.
            n_cores: Number of parallel tasks to run. If <=0, uses all available CPUs.
        """
        self.task_class = task_class
        self.n_cores = (
            n_cores if n_cores > 0 else int(ray.available_resources().get("CPU", 1))
        )
        self.futures: list[ray.ObjectRef] = []
        self.results: list[Any] = []
        self.save_func: Callable[[Any, ...], None] | None = None
        self.save_interval: int = 1

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
        save_func: Callable[[Any, ...], None] | None = None,
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


class RayTaskGeneral:
    """A parent class of others that governs what calculations are run on each
    task."""

    def __init__(self):
        self.results: list[Any] = []

    def run(self, item):
        """Processes a single item.

        Args:
            item: The input data required for the calculation.

        Returns:
            The result of the calculation.
        """
        try:
            result = self.value_func(item)
            return result
        except Exception as e:
            logger.exception(f"Error processing item {item}: {e}")
            return ("error", str(e))

    def value_func(self, item):
        """The definition that actually does the work."""
        pass
