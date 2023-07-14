from complex import *


if __name__ == "__main__":

    dq1 = DualQuaternion(
        Quaternion(1.0, Vector(1.1, 1.2, 0.5)),
        Quaternion(0.5, Vector(0.2, 0.1, 0.3))
    )

    dq2 = DualQuaternion(
        Quaternion(0.1, Vector(0.1, 0.2, 1.5)),
        Quaternion(0.2, Vector(1.2, 1.1, 1.3))
    )

    dq1 *= dq2

    print(dq1)
