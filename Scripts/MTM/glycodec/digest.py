from glypy import plot
from glypy.enzyme import make_n_glycan_pathway

from matplotlib.patches import FancyArrow
from matplotlib import pyplot as plt

glycosylases, glycosyltransferases, seeds = make_n_glycan_pathway()
man9glc3 = seeds[0]


parent_dt, ax = plot.plot(man9glc3)
parent_dt.transform.rotate_deg(-45)

# ----------------------

enz = glycosylases['glucosidaseII']
products = next(enz(man9glc3))

product_dt, ax = plot.plot(products[0], ax=ax)
product_dt.transform.rotate_deg(-45).translate(8, 0)

arrow = FancyArrow(-0.2, 3, 3.8, 0, head_width=1, width=0.5,
                   overhang=0.2, head_length=1, fc=(0.0, 0.5, 0.6))
ax.add_patch(arrow)
ax.text(1., 2.875, "EC:%s" % enz.identifying_information.ec_number, color='white')

# ----------------------

enz = glycosylases['manI']
products = next(enz(products[0]))

product_dt, ax = plot.plot(products[0], ax=ax)
product_dt.transform.rotate_deg(-45).translate(14, -9)

arrow = FancyArrow(8.3, 0.0, 2.5, -3, head_width=1, width=0.5,
                   overhang=0.2, head_length=1, fc=(0.0, 0.5, 0.6))
ax.add_patch(arrow)
ax.text(8.7, -.8, "EC:%s" % enz.identifying_information.ec_number,
        color='white', rotation=-43., )


# ----------------------

enz = glycosylases['manI']
products = next(enz(man9glc3))

product_dt, ax = plot.plot(products[0], ax=ax)
product_dt.transform.rotate_deg(-45).translate(6, -9)

arrow = FancyArrow(-0.3, 0.0, 2.5, -3, head_width=1, width=0.5, overhang=0.2,
                   head_length=1, fc=(0.0, 0.5, 0.6))
ax.add_patch(arrow)
ax.text(0.1, -.8, "EC:%s" % enz.identifying_information.ec_number,
        color='white', rotation=-43., )

# ---------------------
# implicitly the same as the previous

enz = glycosylases['glucosidaseII']
products = next(enz(products[0]))
product_dt, ax = plot.plot(products[0], ax=ax)
product_dt.transform.rotate_deg(-45).translate(14, -9)

arrow = FancyArrow(5.5, -6, 3.8, 0, head_width=1, width=0.5, overhang=0.2,
                   head_length=1, fc=(0.0, 0.5, 0.6))
ax.add_patch(arrow)
ax.text(6.7, -6.15, "EC:%s" % enz.identifying_information.ec_number, color='white')


# ----------------------

enz = glycosylases['glucosidaseII']
ax.text(-8, -4, "EC:%s = %s" % (enz.identifying_information.ec_number,
        enz.identifying_information.name))

enz = glycosylases['manI']
lines = enz.identifying_information.name.split(' ')
lines = ["{0:>{1}}".format(l, 50)
         if i > 0 else l for i, l in enumerate(lines)]
ax.text(-8, -5, "EC:%s = %s" % (enz.identifying_information.ec_number, '\n'.join(lines)))
# ----------------------

ax.set_xlim(-8.1, 16)
ax.set_ylim(-9, 6)
fig = ax.figure
fig.set_figwidth(12)
fig.set_figheight(6)
plt.savefig("enzyme_graph.pdf", bbox_inches='tight')
