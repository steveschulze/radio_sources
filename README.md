# Are there any known radio sources close to source position?

Tool queries FIRST (1.4 GHz), NVSS (1.4 GHz), VLASS (3 GHz) and ASKAP RACS (887.5 GHz) source catalogues and visualises the results.

**Requires**

```
astropy, astroquery, ligo.skymap and plotsettings
```

**How to use it?**

Modify the object properties in `targets.ascii`, search radius in `Cell[4]` and run the notebook.

# Check visibility?

**Requires**

```
astroplan and plotsettings
```

**How to use it?**

Modify the object properties in `targets.ascii` and specify the date in `Cell[2]`. Run the notebook.