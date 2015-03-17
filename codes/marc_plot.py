import matplotlib.pyplot as plt
import numpy as np


def plot(wind1='',wind2='',title1=r'Avg[circle mark] vs. Bl-sinc-W$100\times 50$ [triangle mark], FoV=2$^\circ$',title2=r'Avg[circle mark] vs. Bl-$J_0$-W$1\times 50$',xlabel=r'($\Delta t$[s] ,$\Delta \nu$[Khz])',ylabel='Flux [Jy]'):

# Simple data to display in various forms
    x = np.linspace(0, 2 * np.pi, 400)
    y = np.sin(x ** 2)
    #f, axarr = plt.subplots(2, 2)
    
   # plt.subplot(2, 1, 1)
    StepFreq=np.arange(125000,810000,125000/2.)/1000
    Lav = np.load(wind1+"/Lav0.091667.npy")
    Lmeter = np.load(wind1+"/Lmeter0.091667.npy")
    Rav = (np.load(wind1+"/Rmeter1.750000.npy")*100)#*50)/1000000
    plt.plot(Rav, Lav,'b-o')#,label=r'$r=0,09^{\circ}$')
    plt.plot(Rav, Lmeter,'b-^')
    Lav = np.load(wind1+"/Lav0.250000.npy")
    Lmeter = np.load(wind1+"/Lmeter0.250000.npy")
    plt.plot(Rav, Lav,'y-o')#,label=r'$r=0,25^{\circ}$')
    plt.plot(Rav, Lmeter,'y-^')
    Lav = np.load(wind1+"/Lav0.500000.npy")
    Lmeter = np.load(wind1+"/Lmeter0.500000.npy")
    plt.plot(Rav, Lav,'g-o')#,label=r'$r=0,25^{\circ}$')
    plt.plot(Rav, Lmeter,'g-^')
    Lav = np.load(wind1+"/Lav1.000000.npy")
    Lmeter = np.load(wind1+"/Lmeter1.000000.npy")
    plt.plot(Rav, Lav,'c-o')#,label=r'$r=1^{\circ}$ ')
    plt.plot(Rav, Lmeter,'c-^')
    Lav = np.load(wind1+"/Lav1.500000.npy")
    Lmeter = np.load(wind1+"/Lmeter1.500000.npy")
    plt.plot(Rav, Lav,'k-o')#,label=r'$r=1,5^{\circ}$')
    plt.plot(Rav, Lmeter,'k-^')
    Lav = np.load(wind1+"/Lav2.000000.npy")
    Lmeter = np.load(wind1+"/Lmeter2.000000.npy")
    plt.plot(Rav, Lav,'r-o')#,label=r'$r=2^{\circ}$')
    plt.plot(Rav, Lmeter,'r-^')
  #  Lav = np.load(wind1+"/Lav3.000000.npy")
  #  Lmeter = np.load(wind1+"/Lmeter3.000000.npy")
 #   plt.plot(Rav, Lav,'y-o')#,label=r'$r=2^{\circ}$')
  #  plt.plot(Rav, Lmeter,'y-^')
    #plt.xlim(0, 600)#600)#65)
    plt.ylim(0,1.2)
    plt.grid()
    x = [0]
    y = [0]
    plt.plot(x,y,'b',label=r'$r=0.09^{\circ}$')
    plt.plot(x,y,'y',label=r'$r=0.25^{\circ}$')
    plt.plot(x,y,'g',label=r'$r=0.5^{\circ}$')
    plt.plot(x,y,'c',label=r'$r=1^{\circ}$')
    plt.plot(x,y,'k',label=r'$r=1.5^{\circ}$')
    plt.plot(x,y,'r',label=r'$r=2^{\circ}$')
    
    plt.title(title1,fontsize=18)
    Rav1 = ["" for i in Rav] 
    for i in range(len(Rav)):
         if i%2:
            Rav1[i]=""
         else:Rav1[i]='(%.0f,%.0f)'%(Rav[i],StepFreq[i])       #Rav1[i]="%.1f"%Rav[i] 
			#Rav1[i]='(%.0f,%.0f)'%(Rav[i],StepFreq[i])
    y = np.arange(0,1.01,0.1)
    plt.yticks(y)
    plt.xticks(Rav,Rav1)
    #plt.xticks(Rav,Rav1,visible=False)
    plt.ylabel(ylabel,fontsize=20)
   # plt.legend(loc='down center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=6)
    plt.legend(prop={'size':14},loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=3, fancybox=True, shadow=False)
   # plt.legend(prop={'size':14},loc='upper center', bbox_to_anchor=(.91, 1.02))
    plt.xlabel(xlabel,fontsize=20)
    # plt.subplot(2, 1, 2)
    # Lav = np.load(wind2+"/Lav0.091667.npy")
    # Lmeter = np.load(wind2+"/Lmeter0.091667.npy")
    # Rav = np.load(wind2+"/Rmeter1.750000.npy")/1000
    # plt.plot(Rav, Lav,'b-o',label='r=0,09 deg')
    # plt.plot(Rav, Lmeter,'b-^')
    # Lav = np.load(wind2+"/Lav0.250000.npy")
    # Lmeter = np.load(wind2+"/Lmeter0.250000.npy")
    # plt.plot(Rav, Lav,'g-o',label='r=0,25 deg')
    # plt.plot(Rav, Lmeter,'g-^')
    # Lav = np.load(wind2+"/Lav1.000000.npy")
    # Lmeter = np.load(wind2+"/Lmeter1.000000.npy")
    # plt.plot(Rav, Lav,'c-o',label='r=1 deg')
    # plt.plot(Rav, Lmeter,'c-^')
    # Lav = np.load(wind2+"/Lav1.500000.npy")
    # Lmeter = np.load(wind2+"/Lmeter1.500000.npy")
    # plt.plot(Rav, Lav,'k-o',label='r=1,5 deg')
    # plt.plot(Rav, Lmeter,'k-^')
    # Lav = np.load(wind2+"/Lav2.000000.npy")
    # Lmeter = np.load(wind2+"/Lmeter2.000000.npy")
    # plt.plot(Rav, Lav,'r-o',label='r=2 deg')
    # plt.plot(Rav, Lmeter,'r-^')
    # Lav = np.load(wind2+"/Lav3.000000.npy")
    # Lmeter = np.load(wind2+"/Lmeter3.000000.npy")
    # plt.plot(Rav, Lav,'y-o',label='r=2 deg')
    # plt.plot(Rav, Lmeter,'y-^')
    # plt.xlim(0,1300)
    # plt.ylim(0,1.2)
    # plt.title(title2)
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    # for i in range(len(Rav)):
    #     if i%2:
    #         Rav1[i]=""
    #     else: Rav1[i]=int(Rav[i])
    # plt.xticks(Rav,Rav1)
    # plt.grid()
 #    axarr[0, 0].set_title('Axis [0,0]')
#     axarr[0, 1].plot(Rav, Lav,'b-o',label='r=2 deg')
#     axarr[0, 1].plot(Rav, Lmeter,'b-^')
#     axarr[0, 1].set_title('Axis [0,1]')
#     axarr[1, 0].plot(x, y ** 2,label='here')
#     axarr[1, 0].set_title('Axis [1,0]')
#     axarr[1, 1].plot(x, y ** 2,label='here')
#     axarr[1, 1].set_title('Axis [1,1]')
# # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
    # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
# #box = axarr.get_position()
# #axarr[0,0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
# #axarr[0,1].legend(loc='upper center', bbox_to_anchor=(1, 1.1),
#           fancybox=True, shadow=True, nlig=5)
   
    plt.show()
