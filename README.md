# Multicompartmental-NeuroML
NeuroML branch of the multicompartmental IO model found on: https://github.com/MRIO/IONML & http://www.opensourcebrain.org/projects/io1 


##Network information

Networks are the constructed cell experiments that are able to run in OSB

Networks contains all cells without Q10 parameters
Networks_Q10of3 contains all cells with Q10 factors set to 3
Networks_variableQ10s contains all cells with Q10 factors in accordance with Naomis exploration on https://github.com/njhulst/IO-temperature-dependence




 ## Channel Information
 
 
 All channels are validated against NeuroML_v2.1

Exponential Q10 temperature factor is present in all channels found in channelsQ10 except:
- Cah: `<q10Settings type="q10Fixed" fixedQ10="0.2"/><!-- hack to get 5*tau -->` **Changing Q10 here breaks Calcium dynamics**
- BK


#

 Parameter space analysis needed for (read: run without):
- [ ] BK
- [x] CaCC
- [ ] Cav3.1
