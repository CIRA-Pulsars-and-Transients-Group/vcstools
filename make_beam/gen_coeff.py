from scipy.signal import firwin,kaiser
from pylab import figure,plot,grid,show,xlabel,ylabel,title

N = 128
D = 12
sample_rate = 1.0/1.28E6

cof_bit = 12
b = firwin(N*D-1,1.0/N)
b = (1-pow(2,-(cof_bit-1)))*b/max(b); 
window = kaiser(N*D-1,5)
x_axis=[]
for x in range(-1*(N*D/2),(N*D/2 -1)):
	x_axis.append(x*sample_rate*1E6)

figure(1)
#plot(x_axis,b)
#plot(x_axis,window,'r-')
filter = b*window

#invert the filter
inv_filter = []
threshold = 100;

for coeff in filter:
	inv_coeff = 1.0/coeff;
	if (abs(inv_coeff) < threshold):
		inv_filter.append(inv_coeff)
	else:
		inv_filter.append(threshold*abs(inv_coeff)/(inv_coeff))
        print coeff
plot(x_axis,filter,'g')
#plot(x_axis,inv_filter,'g')
xlabel('Time (us)')
ylabel('Filter/Window coefficients')
title('Fine PFB filter response')
show()

