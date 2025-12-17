[![OpenFOAM](https://img.shields.io/badge/foam_extend_5.0-yellow.svg)](https://sourceforge.net/p/foam-extend/foam-extend-5.0/ci/master/tree/)
[![License](https://img.shields.io/badge/License-GPL--3.0-red.svg)](LICENSE)


# 💻 FSI-IBM Solvers Based on foam-extend-5.0

🚀 A hybrid implementation of **loosely-coupled** solver for **Fluid-Structure Interaction (FSI)** using the **Immersed Boundary Method (IBM)** — built on top of **foam-extend-5.0**.

---

## 👤 Author

**Author:** Tesbo (aka seeeeeeeeeeer)
**Email:** cotsqa@qq.com
**GitHub:** [https://github.com/Tes-bo](https://github.com/Tes-bo)
**Last Updated:** December 17, 2025

---

## 📦 Overview

This project provides an advanced solver for simulating complex FSI problems with IBM:

- 🔄 **Loosely-Coupled FSI Solver**: A more flexible, faster approach for problems where full convergence at every time step isn't necessary.

The **Immersed Boundary Method (IBM)** allows simulation of complex and moving boundaries without the need for body-fitted meshes — perfect for FSI applications.

---

## 🛠️ Features

- ✅ Built on top of `foam-extend-5.0`
- 🧩 Support for **loose** FSI coupling strategies
- 🌊 Immersed Boundary Method for flexible body representation
- 🔧 Modular structure for easy customization and extension
- 📈 Suitable for a wide range of FSI problems: biomedical, aeroelasticity, engineering design, etc.

---

## 📁 Directory Structure

> TODO

---

## 🚀 Getting Started

### 🔧 Requirements

- `foam-extend-5.0` (properly compiled and sourced)
- C++ compiler (e.g. `gcc`, `clang`)
- CMake (if using external libraries)

⚠️ **Note:** This solver requires modifications to the underlying foam-extend-5.0 source code to work properly.

### 🧪 Compilation

```bash
...
wmake
```

---

🧪 Tutorials

Check the tutorials/ directory for sample test cases:
```
cd ...
./Allrun
```

---
🤝 Contribution

Feel free to open issues or submit pull requests! 🛠️
Bug reports, feature requests, and community contributions are welcome.

---

📬 Contact

For questions, suggestions or collaborations:
📧 **Email:** cotsqa@qq.com

---

📝 License

This project is released under the **GPL License**. See LICENSE file for details.

**Citation and Attribution:**
If you use or reference this work, please acknowledge the repository source and the author. Please comply with the GPL license terms.


