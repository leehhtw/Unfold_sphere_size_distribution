# Unfolding sphere size distribution

Here, we implement the recovery of sphere radius histogram in 3d from a circular cross-section radius histogram in 2d. The forward transformation (from 3d to 2d) is well-explained in Eq. (3) of [Taylor, CC, 1983](https://doi.org/10.1111/j.1365-2818.1983.tb04708.x).

The inverse problem (2d to 3d histogram) is solved by using a simple gradient descent.

## Author
* [Hong-Hsi Lee](http://www.diffusion-mri.com/people/hong-hsi-lee)
