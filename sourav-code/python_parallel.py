import multiprocessing as mp
from multiprocessing import Pool
import numpy as np

print(f"Number of cpus: {mp.cpu_count()}")


def worker(num):
    """worker function"""
    print(f"Worker {num}")
    return

if __name__ == '__main__':
    jobs = []
    for i in range(5):
        p = mp.Process(target=worker, args=(i,))
        jobs.append(p)
        p.start()
