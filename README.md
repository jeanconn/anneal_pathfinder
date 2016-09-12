## Realtime tools for the ACA anneal pathfinder procedure

`make_sim_pixel_values.py`: Create simulated pixel data streamed to a file
`plot_pixel_values.py`: Plot the same pixel values
`splat.pl`: version of splat updated to put out pixel data as table in realtime

Simulation Example
------------------

```
$ make_sim_pixel_values.py --delay=0.2

# In a different window plot the brightest 16 pixels in real time
$ plot_pixel_values.py --n-brightest=16

# Plot all 128 pixels in two windows (kinda slow to start)
$ plot_pixel_values.py --pix-filename pixel_values1.dat
$ plot_pixel_values.py --pix-filename pixel_values2.dat
```

Operational Process
-------------------

```
# In one window, run splat to see the hdr3 info
$ memshar_get.pl | ./splat.pl -noimage

# In another window, run splat to get the data tables and see the image data
$ memshar_get.pl | ./splat.pl -noraw -out=6,7

# And in two others, run fits on the brightest pixels in slots 6 and 7
$ plot_pixel_values.py --pix-filename acaslot6.dat
$ plot_pixel_values.py --pix-filename acaslot7.dat
```

