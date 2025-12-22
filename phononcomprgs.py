import numpy as np
import matplotlib.pyplot as plt

data0 = np.load("phononJ0=10^{-3}bk.npz")
gof10,Qlr0,Qphs0,Id0,Ile0,Ire0,Iphs0,cohes0,concv0,Nls0 = data0["gof1"], data0["Qlr"], data0["Qphs"], data0["Id"], data0["Ile"], data0["Ire"], data0["Iphs"], data0["cohes"], data0["concv"], data0["Nls"]

data0gb = np.load("phononJ0=10^{-3}bkglb.npz")
gof10g,Qlr0g,Qphs0g,Id0g,Ile0g,Ire0g,Iphs0g,cohes0g,concv0g,Nls0g = data0gb["gof1"], data0gb["Qlr"], data0gb["Qphs"], data0gb["Id"], data0gb["Ile"], data0gb["Ire"], data0gb["Iphs"], data0gb["cohes"], data0gb["concv"], data0gb["Nls"]


data0red = np.load("phononJ0=10^{-3}bkred.npz")
gof10r,Qlr0r,Qphs0r,cohes0r,concv0r,Nls0r = data0red["gof1"], data0red["Qlrs"], data0red["Qphs"], data0red["cohev"], data0red["concuv"], data0red["Nls"]



data1 = np.load("phononJ0=10^{-2}bk.npz")
gof11,Qlr1,Qphs1,Id1,Ile1,Ire1,Iphs1,cohes1,concv1,Nls1 = data1["gof1"], data1["Qlr"], data1["Qphs"], data1["Id"], data1["Ile"], data1["Ire"], data1["Iphs"], data1["cohes"], data1["concv"], data1["Nls"]

data1gb = np.load("phononJ0=10^{-2}bkglb.npz")
gof11g,Qlr1g,Qphs1g,Id1g,Ile1g,Ire1g,Iphs1g,cohes1g,concv1g,Nls1g = data1gb["gof1"], data1gb["Qlr"], data1gb["Qphs"], data1gb["Id"], data1gb["Ile"], data1gb["Ire"], data1gb["Iphs"], data1gb["cohes"], data1gb["concv"], data1gb["Nls"]

data1red = np.load("phononJ0=10^{-2}bkred.npz")
gof11r,Qlr1r,Qphs1r,cohes1r,concv1r,Nls1r = data1red["gof1"], data1red["Qlrs"], data1red["Qphs"], data1red["cohev"], data1red["concuv"], data1red["Nls"]


data2 = np.load("phononJ0=10^{-1}bk.npz")
gof12,Qlr2,Qphs2,Id2,Ile2,Ire2,Iphs2,cohes2,concv2,Nls2 = data2["gof1"], data2["Qlr"], data2["Qphs"], data2["Id"], data2["Ile"], data2["Ire"], data2["Iphs"], data2["cohes"], data2["concv"], data2["Nls"]

data2gb = np.load("phononJ0=10^{-1}bkglb.npz")
gof12g,Qlr2g,Qphs2g,Id2g,Ile2g,Ire2g,Iphs2g,cohes2g,concv2g,Nls2g = data2gb["gof1"], data2gb["Qlr"], data2gb["Qphs"], data2gb["Id"], data2gb["Ile"], data2gb["Ire"], data2gb["Iphs"], data2gb["cohes"], data2gb["concv"], data2gb["Nls"]

data2red = np.load("phononJ0=10^{-1}bkred.npz")
gof12r,Qlr2r,Qphs2r,cohes2r,concv2r,Nls2r = data2red["gof1"], data2red["Qlrs"], data2red["Qphs"], data2red["cohev"], data2red["concuv"], data2red["Nls"]


data3 = np.load("phononJ0=bk.npz")
gof13,Qlr3,Qphs3,Id3,Ile3,Ire3,Iphs3,cohes3,concv3,Nls3 = data3["gof1"], data3["Qlr"], data3["Qphs"], data3["Id"], data3["Ile"], data3["Ire"], data3["Iphs"], data3["cohes"], data3["concv"], data3["Nls"]

data3gb = np.load("phononJ0=bkglb.npz")
gof13g,Qlr3g,Qphs3g,Id3g,Ile3g,Ire3g,Iphs3g,cohes3g,concv3g,Nls3g = data3gb["gof1"], data3gb["Qlr"], data3gb["Qphs"], data3gb["Id"], data3gb["Ile"], data3gb["Ire"], data3gb["Iphs"], data3gb["cohes"], data3gb["concv"], data3gb["Nls"]

data3red = np.load("phononJ0=bkred.npz")
gof13r,Qlr3r,Qphs3r,cohes3r,concv3r,Nls3r = data3red["gof1"], data3red["Qlrs"], data3red["Qphs"], data3red["cohev"], data3red["concuv"], data3red["Nls"]


data4 = np.load("phononJ0=10bk.npz")
gof14,Qlr4,Qphs4,Id4,Ile4,Ire4,Iphs4,cohes4,concv4,Nls4 = data4["gof1"], data4["Qlr"], data4["Qphs"], data4["Id"], data4["Ile"], data4["Ire"], data4["Iphs"], data4["cohes"], data4["concv"], data4["Nls"]

data4gb = np.load("phononJ0=10bkglb.npz")
gof14g,Qlr4g,Qphs4g,Id4g,Ile4g,Ire4g,Iphs4g,cohes4g,concv4g,Nls4g = data4gb["gof1"], data4gb["Qlr"], data4gb["Qphs"], data4gb["Id"], data4gb["Ile"], data4gb["Ire"], data4gb["Iphs"], data4gb["cohes"], data4gb["concv"], data4gb["Nls"]

data4red = np.load("phononJ0=10bkred.npz")
gof14r,Qlr4r,Qphs4r,cohes4r,concv4r,Nls4r = data4red["gof1"], data4red["Qlrs"], data4red["Qphs"], data4red["cohev"], data4red["concuv"], data4red["Nls"]



data5 = np.load("phononJ0=5_10^{-1}bk.npz")
gof15,Qlr5,Qphs5,Id5,Ile5,Ire5,Iphs5,cohes5,concv5,Nls5 = data5["gof1"], data5["Qlr"], data5["Qphs"], data5["Id"], data5["Ile"], data5["Ire"], data5["Iphs"], data5["cohes"], data5["concv"], data5["Nls"]

data5gb = np.load("phononJ0=5_10^{-1}bkglb.npz")
gof15g,Qlr5g,Qphs5g,Id5g,Ile5g,Ire5g,Iphs5g,cohes5g,concv5g,Nls5g = data5gb["gof1"], data5gb["Qlr"], data5gb["Qphs"], data5gb["Id"], data5gb["Ile"], data5gb["Ire"], data5gb["Iphs"], data5gb["cohes"], data5gb["concv"], data5gb["Nls"]

data5red = np.load("phononJ0=5_10^{-1}bkred.npz")
gof15r,Qlr5r,Qphs5r,cohes5r,concv5r,Nls5r = data5red["gof1"], data5red["Qlrs"], data5red["Qphs"], data5red["cohev"], data5red["concuv"], data5red["Nls"]


plt.figure(figsize=(10, 6))
plt.plot(gof10, Qlr0, lw = 3,label='J/bk=1e-3 (semi)', color='blue')
plt.plot(gof10g, Qlr0g,lw = 3, label='J/bk=1e-3 (glb)', color='red')
plt.plot(gof10r, Qlr0r, lw = 3,linestyle = '--',label='J/bk=1e-3 (red)', color='black')
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$J_{LR}$',fontsize = 20) 
plt.xscale('log')
plt.show()

plt.plot(gof10, Qlr2, lw = 3,label='J/bk=1e-1 (semi)', color='blue')
plt.plot(gof10g, Qlr2g,lw = 3, label='J/bk=1e-1 (glb)', color='red')
plt.plot(gof10r, Qlr2r, lw = 3,linestyle = '--',label='J/bk=1e-1 (red)', color='black')
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$J_{LR}$',fontsize = 20) 
plt.xscale('log')
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(gof10, Qphs0, lw = 3,label='J/bk=1e-3 (semi)', color='blue')
plt.plot(gof10g, Qphs0g, lw = 3,label='J/bk=1e-3 (glb)', color='red')
plt.plot(gof10r, Qphs0r, lw = 3,linestyle = '--',label='J/bk=1e-3 (red)', color='black')
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$J_{ph}$',fontsize = 20) 
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale('log')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(gof10, Qphs2, lw = 3,label='J/bk=1e-1 (semi)', color='blue')
plt.plot(gof10g, Qphs2g, lw = 3,label='J/bk=1e-1 (glb)', color='red')
plt.plot(gof10r, Qphs2r, lw = 3,linestyle = '--',label='J/bk=1e-1 (red)', color='black')
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$J_{ph}$',fontsize = 20) 
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale('log')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(gof10, cohes0, lw = 3,label='J/bk=1e-3 (semi)', color='blue')
plt.plot(gof10g, cohes0g,lw = 3, label='J/bk=1e-3 (glb)', color='red')
plt.plot(gof10r, cohes0r, lw = 3,linestyle = '--',label='J/bk=1e-3 (red)', color='black')
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 20) 
plt.xscale('log')
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(gof10, cohes2, lw = 3,label='J/bk=1e-1 (semi)', color='blue')
plt.plot(gof10g, cohes2g,lw = 3, label='J/bk=1e-1 (glb)', color='red')
plt.plot(gof10r, cohes2r, lw = 3,linestyle = '--',label='J/bk=1e-1 (red)', color='black')
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 20) 
plt.xscale('log')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(gof10, concv0, lw = 3,label='J/bk=1e-3 (semi)', color='blue')
plt.plot(gof10g, concv0g,lw = 3, label='J/bk=1e-3 (glb)', color='red')
plt.plot(gof10r, concv0r, lw = 3,linestyle = '--',label='J/bk=1e-3 (red)', color='black')
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$\mathcal{C}_{on}$',fontsize = 20) 
plt.xscale('log')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(gof10, concv2, lw = 3,label='J/bk=1e-1 (semi)', color='blue')
plt.plot(gof10g, concv2g,lw = 3, label='J/bk=1e-1 (glb)', color='red')
plt.plot(gof10r, concv2r, lw = 3,linestyle = '--',label='J/bk=1e-1 (red)', color='black')
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20) 
plt.ylabel(r'$\mathcal{C}_{on}$',fontsize = 20) 
plt.xscale('log')
plt.show()
