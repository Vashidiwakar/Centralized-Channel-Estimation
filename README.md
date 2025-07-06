# Centralized Channel Estimation in Massive MIMO systems

This repository provides a MATLAB implementation of centralized channel estimation in massive MIMO systems. The algorithm is inspired by the paper:

> **Rajoriya, A., Budhiraja, R., & Hanzo, L. (2022).**  
> _Centralized and decentralized channel estimation in FDD multi-user massive MIMO systems_,  
> IEEE Transactions on Vehicular Technology, 71(7), 7325–7342.  
> [Read Paper on IEEE Xplore](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9749955)

This work focuses only on the centralized estimation part of the paper and exploits structured sparsity in the angular domain, accounting for user-specific, cluster, and common sparsity among scatterers.

---

## File Structure

| File                  | Description |
|-----------------------|-------------|
| `CE_mMIMO.m`          | Main script that runs the entire simulation. It calls all other functions and controls experiment flow. |
| `Channel_Generation.m`| Generates sparse angular-domain channels considering user-specific, cluster-based, and common sparsity. |
| `signal_gen_mMIMO.m`  | Simulates pilot transmission and received signal generation for massive MIMO systems. |
| `create_clusters.m`   | Introduces cluster sparsity in the angular domain. |
| `updates.m`           | Implements the Expectation-Maximization (EM) algorithm for parameter updates. |

---

## Author

**Vashi Diwakar**  
Third-Year Undergraduate, Department of Electrical Engineering, IIT Kanpur  
LinkedIn: [https://www.linkedin.com/in/vashi-diwakar-1b6237288](https://www.linkedin.com/in/vashi-diwakar-1b6237288)  
GitHub: [https://github.com/Vashidiwakar](https://github.com/Vashidiwakar)

---

## Getting Started

### Requirements
- MATLAB R2021a or later
- No special toolboxes are required

### How to Run
Open MATLAB in the project directory and run:

```matlab
>> CE_mMIMO
