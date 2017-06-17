from IPython.display import HTML

HTML('''<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
<form action="javascript:code_toggle()"><input type="submit" value="Click here to toggle on/off the raw code."></form>''')


import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

from scipy.io import loadmat
import ipywidgets as widgets

from IPython import display
w = widgets.Dropdown(
    options={'Pen': 'Pen', 'Pencil': 'Pencil', 'Pad': 'Pad'},
    description='Item:',
)
display.display(w)

def on_button_clicked(b):
    plt.figure(figsize=(20,5))
    plt.plot(range(0,100), np.random.rand(100,1))
    #display.display(plt.gcf())
    display.clear_output(wait=True)
    plt.text(5,0.8, w.value)
    b.description="You clicked"

button = widgets.Button(description="Click Me!")
display.display(button)
button.on_click(on_button_clicked)

# plt.ion()
# fig = plt.figure(figsize=(20,5))
# fig.show()
# fig.canvas.draw()
# ax = fig.add_subplot(1,1,1)