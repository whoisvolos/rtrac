import numpy as np
import numpy.linalg as la
from collections import namedtuple

class Ray3:
    def __init__(self, p0, u):
        self.p0 = np.array(p0)
        npu = np.array(u)
        self.u = npu / la.norm(npu)

    def intr_plane(self, plane):
        denom = np.dot(plane.n, self.u)
        if (denom >= 0):
            return None
        else:
            s = np.dot(plane.n, plane.v0 - ray.p0) / denom
            return self.p0 + s * self.u

class Plane3:
    def __init__(self, v0, n):
        self.v0 = np.array(v0)
        npn = np.array(n)
        self.n = npn / la.norm(npn)

# (p0, u)
ray = Ray3(p0 = [1, 0, 2], u = [2, 0, -1])

# (v0, n)
plane = Plane3(v0 = [1, 0, 0.5], n = [0, 0, 0.5])

for i in range(0, 100000):
    ray.intr_plane(plane)
