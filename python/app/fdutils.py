from dataclasses import dataclass
import sys

import numpy as np
import numpy.typing as npt

import cocoa


@dataclass(kw_only=True)
class Region:
    left: float
    right: float
    bottom: float
    top: float
    eps: float = 1.0e-6

    def inside(self, pt: tuple[float, float]):
        if (self.left - self.eps < pt[0] < self.right + self.eps) and (
            self.bottom - self.eps < pt[1] < self.top + self.eps
        ):
            return True
        return False

    def to_list(self) -> list[float]:
        return [self.left, self.right, self.bottom, self.top]


@dataclass(kw_only=True)
class Target(Region):
    value: float

    def to_list(self):
        return [self.value] + super().to_list()


def quad_interp(x, y, pts):
    xi = np.array([p[0] for p in pts])
    yi = np.array([p[1] for p in pts])
    vi = np.array([p[2] for p in pts])

    A = np.column_stack(
        [
            xi**2 * yi**2,
            xi**2 * yi,
            xi * yi**2,
            xi**2,
            yi**2,
            xi * yi,
            xi,
            yi,
            np.ones_like(xi),
        ]
    )

    coeffs = np.linalg.solve(A, vi)

    value = (
        coeffs[0] * x**2 * y**2
        + coeffs[1] * x**2 * y
        + coeffs[2] * x * y**2
        + coeffs[3] * x**2
        + coeffs[4] * y**2
        + coeffs[5] * x * y
        + coeffs[6] * x
        + coeffs[7] * y
        + coeffs[8]
    )

    return value


@dataclass(kw_only=True)
class TargetFunc(Region):
    values: npt.NDArray

    def __post_init__(self):
        assert len(self.values) == 9
        self.compute_coeffs()

    def to_list(self):
        return self.values.tolist() + super().to_list()

    def pt(self, k: int):
        match k:
            case 0:
                return (self.left, self.bottom)
            case 1:
                return (self.right, self.bottom)
            case 2:
                return (self.right, self.top)
            case 3:
                return (self.left, self.top)
            case 4:
                return (0.5 * (self.left + self.right), self.bottom)
            case 5:
                return (self.right, 0.5 * (self.bottom + self.top))
            case 6:
                return (0.5 * (self.left + self.right), self.top)
            case 7:
                return (self.left, 0.5 * (self.bottom + self.top))
            case 8:
                return (0.5 * (self.left + self.right), 0.5 * (self.bottom + self.top))
            case _:
                sys.exit(1)

    def compute_coeffs(self):
        data = [
            (*self.pt(0), self.values[0]),
            (*self.pt(1), self.values[1]),
            (*self.pt(2), self.values[2]),
            (*self.pt(3), self.values[3]),
            (*self.pt(4), self.values[4]),
            (*self.pt(5), self.values[5]),
            (*self.pt(6), self.values[6]),
            (*self.pt(7), self.values[7]),
            (*self.pt(8), self.values[8]),
        ]
        xi = np.array([p[0] for p in data])
        yi = np.array([p[1] for p in data])
        vi = np.array([p[2] for p in data])

        A = np.column_stack(
            [
                xi**2 * yi**2,
                xi**2 * yi,
                xi * yi**2,
                xi**2,
                yi**2,
                xi * yi,
                xi,
                yi,
                np.ones_like(xi),
            ]
        )

        self.coeffs = np.linalg.solve(A, vi)

    def interp(self, pt: tuple[float, float]) -> float:
        if not self.inside(pt):
            return 0.0
        x, y = pt
        value = (
            self.coeffs[0] * x**2 * y**2
            + self.coeffs[1] * x**2 * y
            + self.coeffs[2] * x * y**2
            + self.coeffs[3] * x**2
            + self.coeffs[4] * y**2
            + self.coeffs[5] * x * y
            + self.coeffs[6] * x
            + self.coeffs[7] * y
            + self.coeffs[8]
        )
        return value


def test_interp():
    p = cocoa.ProblemFD2D()
    p.name = "fd2d_test_target"

    p.mesh = cocoa.MeshFD2D(start=(0.0, 0.0), end=(1.0, 1.0), n=(40, 40))
    p.n_vars = 0
    p.setup_io("output_fd2d_custom")

    target = TargetFunc(
        left=0.5,
        right=0.9,
        bottom=0.1,
        top=0.9,
        values=np.array([1.0, 0.7, 0.6, 0.9, 0.9, 0.8, 0.8, 1.1, 1.0]),
    )

    target_data = cocoa.VectorFD(p.mesh.n_pts, 0.0)
    for j in range(p.mesh.n[1]):
        for i in range(p.mesh.n[0]):
            pt_id = i + p.mesh.n[0] * j
            pt = p.mesh.pt([i, j])
            value = target.interp(pt)
            target_data.set(pt_id, value)
    target = cocoa.FieldMED()
    target.init("target", p, "on_nodes")
    target.set_values(target_data.data)
    target.init_io(p.output_prefix)
    target.print_vtk(0.0, 0)


if __name__ == "__main__":
    test_interp()
