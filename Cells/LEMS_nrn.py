'''
Neuron simulator export for:

Components:
    na_s (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    null (Type: annotation)
    kdr (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    null (Type: annotation)
    k (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    null (Type: annotation)
    cal (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    null (Type: annotation)
    cah (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    kca (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    null (Type: annotation)
    h (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    ca_conc (Type: fixedFactorConcentrationModel:  restingConc=0.0 (SI concentration) decayConstant=0.013333333333333 (SI time) rho=300000.0 (SI rho_factor))
    cacc (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    na_a (Type: ionChannelHH:  conductance=1.0E-11 (SI conductance))
    leak (Type: ionChannelPassive:  conductance=1.0E-11 (SI conductance))
    C51A (Type: cell)
    iClamp0 (Type: pulseGenerator:  delay=0.5 (SI time) duration=0.2 (SI time) amplitude=1.5E-10 (SI current))
    vClamp0 (Type: voltageClampTriple:  active=0.0 (dimensionless) delay=0.05 (SI time) duration=0.2 (SI time) conditioningVoltage=-0.055 (SI voltage) testingVoltage=-0.055 (SI voltage) returnVoltage=-0.055 (SI voltage) simpleSeriesResistance=1000000.0 (SI resistance))
    net (Type: networkWithTemperature:  temperature=308.15 (SI temperature))
    sim1 (Type: Simulation:  length=0.3 (SI time) step=1.0E-5 (SI time))


    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)
         org.neuroml.export  v1.5.5
         org.neuroml.model   v1.5.5
         jLEMS               v0.9.9.2

'''

import neuron

import time
import datetime
import sys

import hashlib
h = neuron.h
h.load_file("stdlib.hoc")

h.load_file("stdgui.hoc")

h("objref p")
h("p = new PythonObject()")

class NeuronSimulation():

    def __init__(self, tstop, dt, seed=123456789):

        print("\n    Starting simulation in NEURON of %sms generated from NeuroML2 model...\n"%tstop)

        self.setup_start = time.time()
        self.report_file = open('report.ex5.txt','w')
        print('Simulator version:  %s'%h.nrnversion())
        self.report_file.write('# Report of running simulation with %s\n'%h.nrnversion())
        self.report_file.write('Simulator=NEURON\n')
        self.report_file.write('SimulatorVersion=%s\n'%h.nrnversion())

        self.report_file.write('SimulationFile=%s\n'%__file__)
        self.report_file.write('PythonVersion=%s\n'%sys.version.replace('\n',' '))
        print('Python version:     %s'%sys.version.replace('\n',' '))
        self.report_file.write('NeuroMLExportVersion=1.5.5\n')
        self.seed = seed
        self.randoms = []
        self.next_global_id = 0  # Used in Random123 classes for elements using random(), etc. 

        self.next_spiking_input_id = 0  # Used in Random123 classes for elements using random(), etc. 

        '''
        Adding simulation Component(id=sim1 type=Simulation) of network/component: net (Type: networkWithTemperature:  temperature=308.15 (SI temperature))
        
        '''

        # Temperature used for network: 308.15 K
        h.celsius = 308.15 - 273.15

        # ######################   Population: pop
        print("Population pop contains 1 instance(s) of component: C51A of type: cell")

        print("Setting the default initial concentrations for ca (used in C51A) to 3.7152 mM (internal), 3.0 mM (external)")
        h("cai0_ca_ion = 3.7152")
        h("cao0_ca_ion = 3.0")

        h.load_file("C51A.hoc")
        a_pop = []
        h("{ n_pop = 1 }")
        h("objectvar a_pop[n_pop]")
        for i in range(int(h.n_pop)):
            h("a_pop[%i] = new C51A()"%i)
            h("access a_pop[%i].Soma"%i)

            self.next_global_id+=1

        h("{ a_pop[0].position(0.0, 800.0, 0.0) }")

        h("proc initialiseV_pop() { for i = 0, n_pop-1 { a_pop[i].set_initial_v() } }")
        h("objref fih_pop")
        h('{fih_pop = new FInitializeHandler(0, "initialiseV_pop()")}')

        h("proc initialiseIons_pop() { for i = 0, n_pop-1 { a_pop[i].set_initial_ion_properties() } }")
        h("objref fih_ion_pop")
        h('{fih_ion_pop = new FInitializeHandler(1, "initialiseIons_pop()")}')

        print("Processing 2 input lists")

        # ######################   Input List: Clamps
        # Adding single input: Component(id=0 type=input)
        h("objref Clamps_0")
        h("a_pop[0].Soma { Clamps_0 = new iClamp0(0.051572595) } ")

        # ######################   Input List: VClamps
        # Adding single input: Component(id=0 type=input)
        h("objref VClamps_0")
        h("a_pop[0].Soma { VClamps_0 = new vClamp0(0.051572595) } ")

        print("Finished processing 2 input lists")

        trec = h.Vector()
        trec.record(h._ref_t)

        h.tstop = tstop

        h.dt = dt

        h.steps_per_ms = 1/h.dt



        # ######################   File to save: time.dat (time)
        # Column: time
        h(' objectvar v_time ')
        h(' { v_time = new Vector() } ')
        h(' { v_time.record(&t) } ')
        h.v_time.resize((h.tstop * h.steps_per_ms) + 1)

        self.initialized = False

        self.sim_end = -1 # will be overwritten

        setup_end = time.time()
        self.setup_time = setup_end - self.setup_start
        print("Setting up the network to simulate took %f seconds"%(self.setup_time))

    def run(self):

        self.initialized = True
        sim_start = time.time()
        print("Running a simulation of %sms (dt = %sms; seed=%s)" % (h.tstop, h.dt, self.seed))

        try:
            h.run()
        except Exception as e:
            print("Exception running NEURON: %s" % (e))
            quit()


        self.sim_end = time.time()
        self.sim_time = self.sim_end - sim_start
        print("Finished NEURON simulation in %f seconds (%f mins)..."%(self.sim_time, self.sim_time/60.0))

        try:
            self.save_results()
        except Exception as e:
            print("Exception saving results of NEURON simulation: %s" % (e))
            quit()


    def advance(self):

        if not self.initialized:
            h.finitialize()
            self.initialized = True

        h.fadvance()


    ###############################################################################
    # Hash function to use in generation of random value
    # This is copied from NetPyNE: https://github.com/Neurosim-lab/netpyne/blob/master/netpyne/simFuncs.py
    ###############################################################################
    def _id32 (self,obj): 
        return int(hashlib.md5(obj.encode('utf-8')).hexdigest()[0:8],16)  # convert 8 first chars of md5 hash in base 16 to int


    ###############################################################################
    # Initialize the stim randomizer
    # This is copied from NetPyNE: https://github.com/Neurosim-lab/netpyne/blob/master/netpyne/simFuncs.py
    ###############################################################################
    def _init_stim_randomizer(self,rand, stimType, gid, seed): 
        #print("INIT STIM  %s; %s; %s; %s"%(rand, stimType, gid, seed))
        rand.Random123(self._id32(stimType), gid, seed)


    def save_results(self):

        print("Saving results at t=%s..."%h.t)

        if self.sim_end < 0: self.sim_end = time.time()


        # ######################   File to save: time.dat (time)
        py_v_time = [ t/1000 for t in h.v_time.to_python() ]  # Convert to Python list for speed...

        f_time_f2 = open('time.dat', 'w')
        num_points = len(py_v_time)  # Simulation may have been stopped before tstop...

        for i in range(num_points):
            f_time_f2.write('%f'% py_v_time[i])  # Save in SI units...
        f_time_f2.close()
        print("Saved data to: time.dat")

        save_end = time.time()
        save_time = save_end - self.sim_end
        print("Finished saving results in %f seconds"%(save_time))

        self.report_file.write('StartTime=%s\n'%datetime.datetime.fromtimestamp(self.setup_start).strftime('%Y-%m-%d %H:%M:%S'))
        self.report_file.write('SetupTime=%s\n'%self.setup_time)
        self.report_file.write('RealSimulationTime=%s\n'%self.sim_time)
        self.report_file.write('SimulationSaveTime=%s\n'%save_time)
        self.report_file.close()

        print("Saving report of simulation to %s"%('report.ex5.txt'))

        print("Done")

        quit()


if __name__ == '__main__':

    ns = NeuronSimulation(tstop=300, dt=0.01, seed=123456789)

    ns.run()
