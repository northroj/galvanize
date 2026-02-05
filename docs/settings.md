# [settings]

Controls global simulation settings and run-time parameters. This block is not repeatable.

| Field | Type | Description |
|--------|------|-------------|
| `num_particles` | int | Number of Monte Carlo particles to simulate |
| `num_t_steps` | int | Number of discrete time steps (Depreciated) |
| `t_step_size` | float | Time step size (shakes) (Depreciated) |
| `dimensions` | int (default = 3) | Affects physics in spatial dimensions (0D-3D). 1D must be x, 2D must be xy.|
| `csd_step` | float (default = 1) | Percentage of the particle's energy to lose per CSD step (1 = 1%) |
| `csd_model` | string | Which model to use for CSD (only spitzer available right now) |
| `scattering_model` | string | Which model to use for scattering. Default = "csd" (for plasma physics using draco). Options: csd_rutherford, gfp2_rutherford (not yet implemented) |
| `verbosity` | int (default = 1) | Enable print statements: 0 = nothing, 10 = all, 7 = particle track level, 4 = timestep level |

---

### Example

```ini

[settings]
num_particles 1000
timestep_bins 0.0 0.01 99 lin
dimensions 3
csd_step 1
csd_model spitzer

```