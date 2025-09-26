#plotting as a bar graph

import matplotlib.pyplot as plt\


x = [1,2,3,4,5]
y = [100, 200,300,400,500]

plt.bar(x,y, color='teal', label='legend')
plt.title( 'graphing attempt')
plt.xlabel( 'x label')
plt.ylabel( 'ylabel')
plt.legend(facecolor='grey', shadow=True, loc=2, title='bar graph')
plt.figtext(.5,.7, "stuff", color="orange", size=25)
plt.show()

image ='C:\\Users\\cnorris19\\Documents\\Class work\\PC .025 3.png'

img = img.imread(image)
plt.imshow(img, extent=(1,5,3,7))
plt.colorbar()
plt.show()

# to make bar plot horizontal make command plt.barh
# to change color make command plt.bar (x,y, color = )
# to add in title insert new line after plt.bar and label plt.title ('insert title')
# to add x or y label insert new line after the plt.title and label plt.xlabel ('x label')
# to insert label add in label= "title for legend" to the plt.plot and then add in a plt.legend line below plt.ylabel
# to change color of legen in plt.legend type in facecolor='color'
# to add shadow to legen type plt.legend( shadow=True)
# can change location of legend by typing plt.legend(loc=number)
# can add legend title by typing plt.legend(title=' insert title')
# to add text in create plt.figtext change ize type size=number
# can plot images in a graph by typing: import matplotlib.image as img, select image that you want and tyupe image= location
# after that type img = img.imread(image), then plt.imshow(img), then finally plt.show
# can add colorbar by typing plt.colorbar
# can control the size of the image by typing extent=(4 numbers) into the plt.imageshow