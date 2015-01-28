#! /usr/bin/env python3
#
# Descriptive statistics!

from math import floor
from time import time

# Decorator interlude.....

def time_this(func):
  def wrapper(*args):
    start_time = time()
    result = func(*args)
    print("Time to find " + func.__name__ + ": " + str(time()-start_time))
    return result
  return wrapper

def indent(space_count):
  def decorator(func):
    def wrapper(*args):
      print(' ' * space_count, end="")
      return func(*args)
    return wrapper
  return decorator

# Actual statistics.....

def print_descriptive_statistics(numbers):
  print("\nDescriptive Statistics: ")
  for key,value in get_descriptive_statistics(numbers).items():
    print("  -" + key + ": " + str(value))
  print("")

def get_descriptive_statistics(numbers):
  stats = {}
  stats["count"] = len(numbers)
  stats["mean"] = mean(numbers)
  stats["median"] = percentile(numbers)
  stats["97.5%"] = percentile(numbers, 0.975)
  stats["2.5%"] = percentile(numbers, 0.025)
  return stats

def mean(numbers):
  return sum(numbers)/len(numbers)

def percentile(numbers, percentile=.5):
  results = []
  sorted_numbers = sorted(numbers)

  # even, average two numbers
  if len(numbers)*percentile % 1.0 == 0.0:
    results.append(sorted_numbers[int(len(numbers)*percentile)])
    results.append(sorted_numbers[int(len(numbers)*percentile - 1)])
  else:
    results.append(sorted_numbers[int(floor(len(numbers)*percentile))])
  return mean(results)


