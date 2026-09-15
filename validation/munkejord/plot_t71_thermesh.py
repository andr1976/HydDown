import os,csv,numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__));navy,red,amber,slate,grey="#002D40","#D61F39","#E6A740","#82979F","#4C4D4E"
th=np.load(os.path.join(H,"t71therm_out.npz")); lu=np.load(os.path.join(H,"test71_dense_out.npz"))
rows=list(csv.reader(open(os.path.join(H,"exp71_blowdown.csv"))));h=rows[0];a=np.array([[float(x) for x in r] for r in rows[1:]]);c={n:i for i,n in enumerate(h)}
tM=a[:,c["t"]];PM=a[:,c["PT163"]];WM=a[:,c["Weight"]];TbM=a[:,c["TT154"]];TtM=a[:,c["TT114"]]
fig,ax=plt.subplots(2,2,figsize=(13,8.5),dpi=150)
ax[0,0].plot(tM,PM,color=slate,lw=2.2,label="measured PT163");ax[0,0].plot(th["t"],th["P"]/1e5,color=red,ls="--",lw=2,label="thermesh wall")
ax[0,0].plot(lu["t"],lu["P"]/1e5,color=amber,ls=":",lw=1.6,label="lumped wall")
ax[0,0].axhline(5.18,color="k",ls=":",lw=0.6);ax[0,0].set_ylabel("pressure (bar)");ax[0,0].set_title("Pressure");ax[0,0].legend(fontsize=8)
ax[0,1].plot(tM,WM,color=slate,lw=2.2,label="measured weight")
ax[0,1].plot(th["t"],th["m"],color=red,ls="--",lw=2,label="thermesh total");ax[0,1].plot(th["t"],th["ms"],color=navy,ls="-.",lw=1.8,label="thermesh dry ice")
ax[0,1].plot(lu["t"],lu["ms"],color=amber,ls=":",lw=1.8,label="lumped dry ice")
ax[0,1].axhline(8.4,color=grey,ls=":",lw=1.2,label="measured retained 8.4 kg")
ax[0,1].set_ylabel("mass (kg)");ax[0,1].set_title("Inventory + dry ice: thermesh vs lumped");ax[0,1].legend(fontsize=7)
ax[1,0].plot(tM,TtM,color=grey,lw=1.6,label="meas top");ax[1,0].plot(tM,TbM,color=slate,lw=2.0,label="meas bottom")
ax[1,0].plot(th["t"],th["Tg"]-273.15,color=red,ls="--",lw=2,label="thermesh GAS");ax[1,0].plot(th["t"],th["Tl"]-273.15,color=navy,ls="-.",lw=2,label="thermesh LIQ/solid")
ax[1,0].set_ylabel("temperature (C)");ax[1,0].set_title("Gas/liquid split (thermesh)");ax[1,0].legend(fontsize=8)
ax[1,1].plot(th["t"],th["ms"],color=red,ls="--",lw=2,label="thermesh dry ice")
ax[1,1].plot(lu["t"],lu["ms"],color=amber,ls=":",lw=2,label="lumped dry ice")
ax[1,1].axhline(8.4,color=grey,ls=":",lw=1.2,label="measured 8.4 kg")
ax[1,1].set_ylabel("dry ice (kg)");ax[1,1].set_title("Retained dry ice: thermesh vs lumped vs measured");ax[1,1].legend(fontsize=8)
for a_ in ax.flat:a_.set_xlabel("time (s)");a_.set_xlim(0,240);a_.grid(alpha=0.3)
fig.suptitle("Munkejord Test 71 (8 mm, no riser) - thermesh vs lumped wall vs 1 Hz data",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.97]);fig.savefig(os.path.join(H,"test71_thermesh_overlay.pdf"));fig.savefig(os.path.join(H,"test71_thermesh_overlay.png"));print("wrote test71_thermesh_overlay",flush=True)
