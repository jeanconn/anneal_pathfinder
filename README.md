## Realtime tools for the ACA anneal pathfinder procedure

`make_sim_pixel_values.py`: Create simulated pixel data streamed to a file
`plot_pixel_values.py`: Plot the same pixel values

Example
-------

```
$ make_sim_pixel_values.py --delay=0.2

# In a different window
$ plot_pixel_values.py --n-brightest=16
```

To Do
-----

- Add sherpa fitting of a dark current profile that uses `t_ccd`
  to individually fit for `dark_current` and `scale_factor` for
  each bright pixel.
- Report fits at each iteration along with logging other useful
  information.

Sherpa fitting
---------------

This is fitting for `scale` in the modified version of
`mica.archive.aca_dark.dark_cal.dark_temp_scale`:

```
def dark_temp_scale(t_ccd, t_ccd_ref=-19.0, scale=0.70):
    """
    Return the multiplicative scale factor to convert a CCD dark map from
    the actual temperature ``t_ccd`` to the reference temperature ``t_ccd_ref``.

    Based on best global fit for dark current model in `plot_predicted_warmpix.py`.
    Previous value was 0.62 instead of 0.70.  This represents the change in
    dark current for each 4 degC decrease::

      >>> from mica.archive.aca_dark import temp_scalefac
      >>> print temp_scalefac(t_ccd=-15, t_ccd_ref=-19)
      0.7

    :param t_ccd: actual temperature (degC)
    :param t_ccd_ref: reference temperature (degC, default=-19.0)

    :returns: scale factor
    """
    return np.exp(np.log(scale) / 4.0 * (t_ccd - t_ccd_ref))
```
