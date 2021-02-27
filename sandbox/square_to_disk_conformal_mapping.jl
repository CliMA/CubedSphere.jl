using CairoMakie
using CubedSphere

w(z) = (1 - cn(z, 1/2)) / sn(z, 1/2)

CubedSphere.visualize_conformal_mapping(w, x_min=-1.85, y_min=-1.85, x_max=1.85, y_max=1.85, n_lines=51, n_samples=101, filepath="square_to_disk_conformal_mapping.png")
