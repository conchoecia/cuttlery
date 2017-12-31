#!/usr/bin/env python

import pandas as pd
import numpy as np

s = pd.Series(np.random.randn(5))
print(getattr(s, "mean")())
print(s.mean())
