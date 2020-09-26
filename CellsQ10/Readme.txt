 ## Channel Information
 
 
 All channels are validated against NeuroML_v2.1

Exponential Q10 temperature factor is present in all channels except:
- Cah: `<q10Settings type="q10Fixed" fixedQ10="0.2"/><!-- hack to get 5*tau -->` **Changing Q10 here breaks Calcium dynamics**
- BK

Experimentaltemp is set **@Roomtemperature** or taken from literature
Roomtemperature = **20degC**
#
#
#
#
#
#
 Parameter space analysis needed for (read: run without):
- [ ] BK
- [ ] CaCC
- [ ] Cav3.1
