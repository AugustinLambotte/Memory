import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime, date, timedelta
import matplotlib.path as mpath
import cmocean as cmo
import cartopy.crs as ccrs


Ricker = np.array([[-267,-21,-540,-279,np.nan,np.nan,np.nan,np.nan,np.nan,-164,-214,-354],
                  [-129,-381,-379,-487,np.nan,np.nan,np.nan,np.nan,np.nan,-203,-182,-187],
                  [-103,-163,-299,-318,np.nan,np.nan,np.nan,np.nan,np.nan,-215,-400,-231],
                  [-78,-195,-345,-452,np.nan,np.nan,np.nan,np.nan,np.nan,-200,-165,-373],
                  [-160,-425,-429,-354,np.nan,np.nan,np.nan,np.nan,np.nan,-52,-261,-275],
                  [-177,-352,-348,-310,np.nan,np.nan,np.nan,np.nan,np.nan,-129,-151,-307]])
        
M1 = np.array([[-238,-24,-478,-255,np.nan,np.nan,np.nan,np.nan,np.nan,-149,-163,-293],
              [-109,-299,-287,-428,np.nan,np.nan,np.nan,np.nan,np.nan,-207,-157,-125],
              [-80,-122,-254,-254,np.nan,np.nan,np.nan,np.nan,np.nan,-212,-372,-211],
              [-49,-105,-240,-401,np.nan,np.nan,np.nan,np.nan,np.nan,-203,-122,-307],
              [-129,-358,-328,-284,np.nan,np.nan,np.nan,np.nan,np.nan,-72,-215,-243],
              [-129,-272,-255,-264,np.nan,np.nan,np.nan,np.nan,np.nan,-98,-90,-243]])
        
M2 = np.array([[-238,-34,-442,-230,-278,-185,-115,-64,-28,-151,-175,-290],
              [-137,-300,-267,-372,-334,-218,-187,-131,-100,-160,-149,-136],
              [-78,-109,-217,-219,-194,-140,-107,-98,-26,-228,-367,-191],
              [-61,-114,-282,-425,-232,-161,-112,-184,-194,-170,-162,-283],
              [-129,-355,-339,-308,-171,-240,-114,-11,-107,-78,-192,-244],
              [-150,-267,-287,-289,-196,-194,-113,-198,-75,-97,-72,-222]])

np.savetxt('Data/Ricker.txt',Ricker)
np.savetxt('Data/M1.txt',M1)
np.savetxt('Data/M2.txt',M2)