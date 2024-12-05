from typing import Any

import multiprocessing
import platform

from loguru import logger


class MultiprocessingManager:
    """A class for running calculations on multiple processors."""

    def __init__(self, inputs, num_processors, task_class):
        """
        Launches a calculation on multiple processors.

        Args:
            inputs: A list containing all the input required for the calculation.
            num_processors: The number of processors to use. If <=0, uses all available.
            task_class: A class with a method to process inputs.
        """
        self.results = []

        if platform.system().upper().startswith("WIN"):
            multiprocessing.freeze_support()

        if num_processors != 1 and platform.system().upper().startswith("WIN"):
            print(
                "WARNING: Use of multiple processors is not fully supported on Windows. Proceeding with one processor..."
            )
            num_processors = 1

        if num_processors == 1:
            # Single processor execution
            task_instance = task_class()
            self.results = [task_instance.runit(item) for item in inputs]
        else:
            # Multiple processors execution
            if num_processors <= 0:
                num_processors = multiprocessing.cpu_count()

            num_processors = min(num_processors, len(inputs)) if len(inputs) > 0 else 1

            if len(inputs) == 0:
                self.results = []
                return

            try:
                with multiprocessing.Pool(processes=num_processors) as pool:
                    # Initialize separate instances for each process
                    task_instances = [task_class() for _ in range(num_processors)]
                    self.results = pool.map(
                        lambda args: args[0].runit(args[1]), zip(task_instances, inputs)
                    )
            except KeyboardInterrupt:
                print("Process interrupted by user. Terminating pool...")
                pool.terminate()
                pool.join()
                raise
            except Exception as e:
                print(f"An error occurred: {e}. Terminating pool...")
                pool.terminate()
                pool.join()
                raise


class MultiprocessingTaskGeneral:
    """A parent class of others that governs what calculations are run on each
    process."""

    def __init__(self):
        self.results: list[Any] = []

    def runit(self, item):
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


class GeneralTask:
    """A class that determines the specific calculations that will be
    performed when multi-processor support is used. Other, more
    specific classes will inherit this one."""

    results = []

    def runit(self, running, mutex, results_queue, items):
        for item in items:
            self.value_func(item, results_queue)
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

    # this is the function that changes through inheritance
    def value_func(self, item, results_queue):
        print(item)  # here's where you do something
        self.results.append(item)  # here save the results for later compilation
