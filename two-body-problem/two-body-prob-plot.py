import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Obtain the simulation name
parser = argparse.ArgumentParser()
parser.add_argument('--name', type=str, help="Input simulation name", required=True)
name = parser.parse_args().name
print("2-body central force problem simulation")
print(f"Plotting for simulation: {name}")

# Check for data file
pwd = os.path.dirname(os.path.abspath(__file__))
data_file_dir = os.path.join(pwd, name)
data_file_loc = os.path.join(pwd, name, name+"_data.txt")

if os.path.isdir(data_file_dir):
	if os.listdir(data_file_dir): ("Data is read from: %s"%(data_file_loc))
	else: sys.exit(f"Error: {data_file_dir} directory is empty.\nExiting code.")
else:
	sys.exit("Error: {data_file_dir} directory does not exist.\nExiting code.")

x, t, r, th, pr = np.loadtxt(data_file_loc, dtype=np.float64).T

plt.figure()
plt.xlabel("Time")
plt.ylabel(r"$r$")
plt.title(r"$r$ v/s $t$")
plt.plot(t,r,ls="-",color="black")
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_r_vs_t.png"))
plt.close()

plt.figure()
plt.xlabel("Time")
plt.ylabel(r"$\theta$")
plt.title(r"$\theta$ v/s $t$")
plt.plot(t,th,ls="-",color="black")
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_theta_vs_t.png"))
plt.close()

plt.figure()
plt.xlabel("Time")
plt.ylabel(r"$p_r$")
plt.title(r"$p_r$ v/s $t$")
plt.plot(t,pr,ls="-",color="black")
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_pr_vs_t.png"))
plt.close()

plt.figure()
plt.xlabel("Time")
plt.title(r"$r$, $\theta$ and $p_r$ v/s $t$")
plt.plot(t,r,ls="-",color="black",label=r"$r$")
plt.plot(t,th,ls="--",color="black",label=r"$\theta$")
plt.plot(t,pr,ls=":",color="black",label=r"$p_r$")
plt.legend(frameon=False)
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_all_vs_t.png"))
plt.close()

x = r*np.cos(th)
y = r*np.sin(th)
plt.figure()
plt.title("Trajectory")
plt.plot(x,y,ls="-",color="black")
plt.scatter([0],[0],marker="x",color="black")
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_trajectory.png"))
plt.close()

plt.figure()
plt.xlabel("Distance Covered")
plt.ylabel(r"$r$")
plt.title(r"$r$ v/s $x$")
plt.plot(x,r,ls="-",color="black")
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_r_vs_x.png"))
plt.close()

plt.figure()
plt.xlabel("Distance Covered")
plt.ylabel(r"$\theta$")
plt.title(r"$\theta$ v/s $x$")
plt.plot(x,th,ls="-",color="black")
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_theta_vs_x.png"))
plt.close()

plt.figure()
plt.xlabel("Distance Covered")
plt.ylabel(r"$p_r$")
plt.title(r"$p_r$ v/s $x$")
plt.plot(x,pr,ls="-",color="black")
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_pr_vs_x.png"))
plt.close()

plt.figure()
plt.xlabel("Distance Covered")
plt.title(r"$r$, $\theta$ and $p_r$ v/s $x$")
plt.plot(x,r,ls="-",color="black",label=r"$r$")
plt.plot(x,th,ls="--",color="black",label=r"$\theta$")
plt.plot(x,pr,ls=":",color="black",label=r"$p_r$")
plt.legend(frameon=False)
plt.savefig(os.path.join(data_file_dir, "plot_"+name+"_all_vs_x.png"))
plt.close()

print("Plotting done")
