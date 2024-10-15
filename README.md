# DRPM: Degeneracy Resilient Point-to-Plane Error Minimization

A probabilistic method to detect and handle degeneracies in point-to-plane error minimization for LiDAR SLAM. 

Presented in the paper *Probabilistic Degeneracy Detection for Point-to-Plane Error Minimization* ([preprint](https://arxiv.org/abs/2410.10784) | [video](https://www.youtube.com/watch?v=bKnHs_wwnXs)).

## Usage

`src/degeneracy.h` contains a reference implementation of the method proposed in the paper "Probabilistic Degeneracy Detection for Point-to-Plane Error Minimization". See the paper for the details. `src/example.cpp` file exemplifies usage. The implementation can be used in your point cloud registration pipeline by including the `degeneracy.h` file in your project.

## Building the example

See the `Dockerfile` for how to build the example.

To build and run the example using docker, run the following commands:

```bash
docker build -t drpm .
docker run --rm -it drpm
```
Docker must be installed on your machine for this to work.

## Citation

If you use this code in your research, please cite the following paper:
```
@misc{hatleskog2024drpm,
  title={Probabilistic Degeneracy Detection for Point-to-Plane Error Minimization},
  author={Hatleskog, Johan and Alexis, Kostas},
  year={2024},
  note={arXiv:2410.10784}
}
```
