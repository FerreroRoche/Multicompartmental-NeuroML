// Generated code block BEGIN
#define M_PI       3.14159265358979323846
#include <math.h>
typedef float * __restrict__ __attribute__((align_value (32))) Table_F32;
typedef long long * __restrict__ __attribute__((align_value (32))) Table_I64;
static float stepf( float x ){ if( x < 0 ) return 0; else return 1;  }

// Credits to Thomas T. Wang: wang@cup.hp.com
static unsigned long long hash64shift( unsigned long long key ){
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}
static unsigned long long hash_128_to_64( unsigned long long hi, unsigned long long lo ){
	return hash64shift( hash64shift( lo ) ^ hi );
}

static float randof( float x, long long work_item, long long instance, long long step, int invocation_id ){
	// Make a unique stamp for the random number sampled
	// Unique factors: work item, tabular instance, serial number of RNG invocation in kernel, timestep 
	// Capacities: 1T work items, 16M instances, 64K invocations, 1T timesteps 
	unsigned long long stamp_hi = work_item * (1ULL << 24) | instance % (1ULL << 24);
	unsigned long long stamp_lo = invocation_id * (1ULL << 40) | step % (1ULL << 40);
	unsigned long long sample = hash_128_to_64( stamp_hi, stamp_lo );
	const/*ant*/int sample_scale = (1 << 23);
	float result = ( (float) ( sample % sample_scale ) ) / ( (float) (sample_scale) );
	return x * result;
}

void doit( double time, float dt, const float *__restrict__ global_constants, long long const_local_index, 
const long long *__restrict__ global_const_table_f32_sizes, const Table_F32 *__restrict__ global_const_table_f32_arrays, long long table_cf32_local_index,
const long long *__restrict__ global_const_table_i64_sizes, const Table_I64 *__restrict__ global_const_table_i64_arrays, long long table_ci64_local_index,
const long long *__restrict__ global_state_table_f32_sizes, const Table_F32 *__restrict__ global_state_table_f32_arrays, Table_F32 *__restrict__ global_stateNext_table_f32_arrays, long long table_sf32_local_index,
const long long *__restrict__ global_state_table_i64_sizes,       Table_I64 *__restrict__ global_state_table_i64_arrays, Table_I64 *__restrict__ global_stateNext_table_i64_arrays, long long table_si64_local_index,
const float *__restrict__ global_state, float *__restrict__ global_stateNext, long long state_local_index, 
long long step ){
	
	
	char initial_state = (step <= 0);
	const float time_f32 = time; //when not accumulating small deltas, double precision is not necessary, and it messes up with SIMD
	
	const long long NOT_AN_INSTANCE = ~0xFee1600dLL; // if it's misused to index an array it will probably stop right there ㋡
	long long instance = NOT_AN_INSTANCE; // for RNG use
	long long rng_offset = 0; // for RNG use too
	
	const long long const_cell_index = const_local_index;
	const long long state_cell_index = state_local_index;
	const long long table_cf32_cell_index = table_cf32_local_index;
	const long long table_ci64_cell_index = table_ci64_local_index;
	const long long table_sf32_cell_index = table_sf32_local_index;
	const long long table_si64_cell_index = table_si64_local_index;
	
	const float *cell_constants = global_constants + const_cell_index;
	const float *cell_state     = global_state     + state_cell_index;
	      float *cell_stateNext = global_stateNext + state_cell_index;
	
		const long long *cell_const_table_f32_sizes      = global_const_table_f32_sizes      + table_cf32_cell_index;
		const Table_F32 *cell_const_table_f32_arrays     = global_const_table_f32_arrays     + table_cf32_cell_index;
		const long long *cell_const_table_i64_sizes      = global_const_table_i64_sizes      + table_ci64_cell_index;
		const Table_I64 *cell_const_table_i64_arrays     = global_const_table_i64_arrays     + table_ci64_cell_index;
		const long long *cell_state_table_f32_sizes      = global_state_table_f32_sizes      + table_sf32_cell_index;
		const Table_F32 *cell_state_table_f32_arrays     = global_state_table_f32_arrays     + table_sf32_cell_index;
		      Table_F32 *cell_stateNext_table_f32_arrays = global_stateNext_table_f32_arrays + table_sf32_cell_index;
		const long long *cell_state_table_i64_sizes      = global_state_table_i64_sizes      + table_si64_cell_index;
		      Table_I64 *cell_state_table_i64_arrays     = global_state_table_i64_arrays     + table_si64_cell_index;
		      Table_I64 *cell_stateNext_table_i64_arrays = global_stateNext_table_i64_arrays + table_si64_cell_index;
	
	const float temperature = cell_constants[1244]; //a global if there ever was one
	
	const float *V = &cell_state[0]; 
	      float *V_next = &cell_stateNext[0]; 
	const float *R_Axial = &cell_constants[311]; 
	const float *C = &cell_constants[0]; 
	const float *V_threshold = &cell_constants[622]; 
	const float *Area = &cell_constants[933]; 
	
	const Table_I64 Comp_Coff    = cell_const_table_i64_arrays[0];
	const Table_I64 Comp_Soff    = cell_const_table_i64_arrays[1];
	const Table_I64 Comp_CF32off = cell_const_table_i64_arrays[2];
	const Table_I64 Comp_SF32off = cell_const_table_i64_arrays[3];
	const Table_I64 Comp_CI64off = cell_const_table_i64_arrays[4];
	const Table_I64 Comp_SI64off = cell_const_table_i64_arrays[5];
	const Table_I64 Comp_Roff    = cell_const_table_i64_arrays[6];
	// Internal Code for compartment type 0
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[7];
	const long long Type_Compartments    = cell_const_table_i64_sizes [7];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	float Acomp = Area[comp];
	float Vcomp = V[comp];
	float I_internal = 0;
	// Ion flux sources
	// Ion concentrations
	const float Ca_concentration = 0;
	const float Ca_concentration_extra = 0;
	const float Ca2_concentration = 0;
	const float Ca2_concentration_extra = 0;
	// Current from ion channels
	float I_channels_total = 0;
	{
	float Vshift = 0;
		float Erev  = local_constants[0];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
		{
		float x; // define exposure
		x = local_constants[1] / (1 + exp( (local_constants[2] - Vcomp ) / local_constants[3] ) );
		chan_gate_0_q = x;
		}
	float chan_gate_1_q; 
	chan_gate_1_q = local_state[0]; 
	// dynamics for channel 0 gate 1 
	{
		float q10 = local_constants[4];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 1 for Fixed channel 0
		float Lems_property_0 = local_constants[5];
		float Lems_property_1 = local_constants[6];
		float Lems_property_2 = local_constants[7];
		// state variables HHRate BaseTau 1 for Fixed channel 0
		// declare derived variables HHRate BaseTau 1 for Fixed channel 0
		float Lems_derived_0 = NAN;
		// common read-only namespace? HHRate BaseTau 1 for Fixed channel 0
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_property_2;
		float *Lems_assigned_4 = &Lems_derived_0;
		// compute derived HHRate BaseTau 1 for Fixed channel 0
		Lems_derived_0 = ( ( *Lems_assigned_1/* time */ * ( expf( ( ( ( ( *Lems_assigned_0/* voltage */ - *Lems_assigned_2/* voltage */ ) )/* voltage */ / *Lems_assigned_3/* voltage */ ) )/* unitless */ ) )/* unitless */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 1 for Fixed channel 0
		float Lems_exposure_t = Lems_derived_0;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[8] / (1 + exp( (local_constants[9] - Vcomp ) / local_constants[10] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[0] = inf;
		}else{
			local_stateNext[0] = local_state[0] + dt * ( ( inf - local_state[0] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_1_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[11]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[12];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[1]; 
	// dynamics for channel 1 gate 0 
	{
		float q10 = local_constants[13];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 0 for Fixed channel 1
		float Lems_property_0 = local_constants[14];
		float Lems_property_1 = local_constants[15];
		// state variables HHRate BaseTau 0 for Fixed channel 1
		// declare derived variables HHRate BaseTau 0 for Fixed channel 1
		float Lems_derived_0 = NAN;
		float Lems_derived_1 = NAN;
		// common read-only namespace? HHRate BaseTau 0 for Fixed channel 1
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_derived_0;
		float *Lems_assigned_4 = &Lems_derived_1;
		// compute derived HHRate BaseTau 0 for Fixed channel 1
		Lems_derived_0 = ( ( *Lems_assigned_0/* voltage */ / *Lems_assigned_2/* voltage */ ) )/* unitless */;
		Lems_derived_1 = ( ( ( ( 5/* unitless */ + ( ( 47/* unitless */ * ( expf( ( ( ( ( - ( ( ( ( - 50/* unitless */ ) )/* unitless */ - *Lems_assigned_3/* unitless */ ) )/* unitless */ ) )/* unitless */ / 900/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ * *Lems_assigned_1/* time */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 0 for Fixed channel 1
		float Lems_exposure_t = Lems_derived_1;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[16] / (1 + exp( (local_constants[17] - Vcomp ) / local_constants[18] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[1] = inf;
		}else{
			local_stateNext[1] = local_state[1] + dt * ( ( inf - local_state[1] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[19]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[20];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[2]; 
	// dynamics for channel 2 gate 0 
	{
		float q10 = local_constants[21];
		float alpha;
		{
		float r; // define exposure
		r = local_constants[22] * ( ( Vcomp == local_constants[23]) ? 1 : ( ( (Vcomp - local_constants[23] ) / local_constants[24] )  / (1 - exp( - (Vcomp - local_constants[23] ) / local_constants[24] ) ) ) );
		alpha = r;
		}
		float beta;
		{
		float r; // define exposure
		r = local_constants[25] * exp( (Vcomp - local_constants[26] ) / local_constants[27] );
		beta = r;
		}
		float tau;
		tau = 1 / ( alpha + beta );
		float inf;
		inf = alpha / ( alpha + beta );
		if(initial_state){
			local_stateNext[2] = inf;
		}else{
			local_stateNext[2] = local_state[2] + dt * ( ( inf - local_state[2] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[28]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[29];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[3]; 
	// dynamics for channel 3 gate 0 
	{
		float q10 = local_constants[30];
		float tau;
		{
		float t; // define exposure
		t = local_constants[31];
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[32] / (1 + exp( (local_constants[33] - Vcomp ) / local_constants[34] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[3] = inf;
		}else{
			local_stateNext[3] = local_state[3] + dt * ( ( inf - local_state[3] ) / tau ) * q10 ;
		}
	}
	float chan_gate_1_q; 
	chan_gate_1_q = local_state[4]; 
	// dynamics for channel 3 gate 1 
	{
		float q10 = local_constants[35];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 1 for Fixed channel 3
		float Lems_property_0 = local_constants[36];
		float Lems_property_1 = local_constants[37];
		// state variables HHRate BaseTau 1 for Fixed channel 3
		// declare derived variables HHRate BaseTau 1 for Fixed channel 3
		float Lems_derived_0 = NAN;
		float Lems_derived_1 = NAN;
		// common read-only namespace? HHRate BaseTau 1 for Fixed channel 3
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_derived_0;
		float *Lems_assigned_4 = &Lems_derived_1;
		// compute derived HHRate BaseTau 1 for Fixed channel 3
		Lems_derived_0 = ( ( *Lems_assigned_0/* voltage */ / *Lems_assigned_2/* voltage */ ) )/* unitless */;
		Lems_derived_1 = ( ( *Lems_assigned_1/* time */ * ( ( ( ( ( ( 20/* unitless */ * ( expf( ( ( ( ( *Lems_assigned_3/* unitless */ + 160/* unitless */ ) )/* unitless */ / 30/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ / ( ( 1/* unitless */ + ( expf( ( ( ( ( *Lems_assigned_3/* unitless */ + 84/* unitless */ ) )/* unitless */ / 7.2999999999999998/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ + 35/* unitless */ ) )/* unitless */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 1 for Fixed channel 3
		float Lems_exposure_t = Lems_derived_1;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[38] / (1 + exp( (local_constants[39] - Vcomp ) / local_constants[40] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[4] = inf;
		}else{
			local_stateNext[4] = local_state[4] + dt * ( ( inf - local_state[4] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_1_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[41]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[42];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling ;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[43]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	// Current from synapses
	float I_synapses_total = 0;
	// Current from inputs
	float I_input_total = 0;
	{
	const float     *Weight     = local_const_table_f32_arrays[0];
	// Pulse inputs
	const long long Instances_input_pulse = local_const_table_f32_sizes[1]; //same for all parallel arrays
	const float     *Imax_input_pulse     = local_const_table_f32_arrays[1];
	const float     *Start_input_pulse    = local_const_table_f32_arrays[2];
	const float     *Duration_input_pulse = local_const_table_f32_arrays[3];
	float I_input_pulse = 0;
	for(long long instance = 0; instance < Instances_input_pulse; instance++){
		if( Start_input_pulse[instance] <= time && time <=  Start_input_pulse[instance] +  Duration_input_pulse[instance] ) I_input_pulse += Imax_input_pulse[instance] * Weight[instance];
	}
	I_input_total += I_input_pulse;

	}
	I_internal = I_channels_total + I_input_total + I_synapses_total;
	if(initial_state){
		// initialize
		V_next[comp] = V[comp];
	}else{
		V_next[comp] = V[comp] + ( dt * ( I_internal ) / C[comp] );
	}
	}
	}
	// Internal Code for compartment type 0 end
	// Internal Code for compartment type 1
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[8];
	const long long Type_Compartments    = cell_const_table_i64_sizes [8];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	float Acomp = Area[comp];
	float Vcomp = V[comp];
	float I_internal = 0;
	// Ion flux sources
	// Ion concentrations
	const float Ca_concentration = 0;
	const float Ca_concentration_extra = 0;
	const float Ca2_concentration = 0;
	const float Ca2_concentration_extra = 0;
	// Current from ion channels
	float I_channels_total = 0;
	{
	float Vshift = 0;
		float Erev  = local_constants[0];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
		{
		float x; // define exposure
		x = local_constants[1] / (1 + exp( (local_constants[2] - Vcomp ) / local_constants[3] ) );
		chan_gate_0_q = x;
		}
	float chan_gate_1_q; 
	chan_gate_1_q = local_state[0]; 
	// dynamics for channel 0 gate 1 
	{
		float q10 = local_constants[4];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 1 for Fixed channel 0
		float Lems_property_0 = local_constants[5];
		float Lems_property_1 = local_constants[6];
		float Lems_property_2 = local_constants[7];
		// state variables HHRate BaseTau 1 for Fixed channel 0
		// declare derived variables HHRate BaseTau 1 for Fixed channel 0
		float Lems_derived_0 = NAN;
		// common read-only namespace? HHRate BaseTau 1 for Fixed channel 0
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_property_2;
		float *Lems_assigned_4 = &Lems_derived_0;
		// compute derived HHRate BaseTau 1 for Fixed channel 0
		Lems_derived_0 = ( ( *Lems_assigned_1/* time */ * ( expf( ( ( ( ( *Lems_assigned_0/* voltage */ - *Lems_assigned_2/* voltage */ ) )/* voltage */ / *Lems_assigned_3/* voltage */ ) )/* unitless */ ) )/* unitless */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 1 for Fixed channel 0
		float Lems_exposure_t = Lems_derived_0;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[8] / (1 + exp( (local_constants[9] - Vcomp ) / local_constants[10] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[0] = inf;
		}else{
			local_stateNext[0] = local_state[0] + dt * ( ( inf - local_state[0] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_1_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[11]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[12];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[1]; 
	// dynamics for channel 1 gate 0 
	{
		float q10 = local_constants[13];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 0 for Fixed channel 1
		float Lems_property_0 = local_constants[14];
		float Lems_property_1 = local_constants[15];
		// state variables HHRate BaseTau 0 for Fixed channel 1
		// declare derived variables HHRate BaseTau 0 for Fixed channel 1
		float Lems_derived_0 = NAN;
		float Lems_derived_1 = NAN;
		// common read-only namespace? HHRate BaseTau 0 for Fixed channel 1
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_derived_0;
		float *Lems_assigned_4 = &Lems_derived_1;
		// compute derived HHRate BaseTau 0 for Fixed channel 1
		Lems_derived_0 = ( ( *Lems_assigned_0/* voltage */ / *Lems_assigned_2/* voltage */ ) )/* unitless */;
		Lems_derived_1 = ( ( ( ( 5/* unitless */ + ( ( 47/* unitless */ * ( expf( ( ( ( ( - ( ( ( ( - 50/* unitless */ ) )/* unitless */ - *Lems_assigned_3/* unitless */ ) )/* unitless */ ) )/* unitless */ / 900/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ * *Lems_assigned_1/* time */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 0 for Fixed channel 1
		float Lems_exposure_t = Lems_derived_1;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[16] / (1 + exp( (local_constants[17] - Vcomp ) / local_constants[18] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[1] = inf;
		}else{
			local_stateNext[1] = local_state[1] + dt * ( ( inf - local_state[1] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[19]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[20];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[2]; 
	// dynamics for channel 2 gate 0 
	{
		float q10 = local_constants[21];
		float alpha;
		{
		float r; // define exposure
		r = local_constants[22] * ( ( Vcomp == local_constants[23]) ? 1 : ( ( (Vcomp - local_constants[23] ) / local_constants[24] )  / (1 - exp( - (Vcomp - local_constants[23] ) / local_constants[24] ) ) ) );
		alpha = r;
		}
		float beta;
		{
		float r; // define exposure
		r = local_constants[25] * exp( (Vcomp - local_constants[26] ) / local_constants[27] );
		beta = r;
		}
		float tau;
		tau = 1 / ( alpha + beta );
		float inf;
		inf = alpha / ( alpha + beta );
		if(initial_state){
			local_stateNext[2] = inf;
		}else{
			local_stateNext[2] = local_state[2] + dt * ( ( inf - local_state[2] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[28]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[29];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[3]; 
	// dynamics for channel 3 gate 0 
	{
		float q10 = local_constants[30];
		float tau;
		{
		float t; // define exposure
		t = local_constants[31];
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[32] / (1 + exp( (local_constants[33] - Vcomp ) / local_constants[34] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[3] = inf;
		}else{
			local_stateNext[3] = local_state[3] + dt * ( ( inf - local_state[3] ) / tau ) * q10 ;
		}
	}
	float chan_gate_1_q; 
	chan_gate_1_q = local_state[4]; 
	// dynamics for channel 3 gate 1 
	{
		float q10 = local_constants[35];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 1 for Fixed channel 3
		float Lems_property_0 = local_constants[36];
		float Lems_property_1 = local_constants[37];
		// state variables HHRate BaseTau 1 for Fixed channel 3
		// declare derived variables HHRate BaseTau 1 for Fixed channel 3
		float Lems_derived_0 = NAN;
		float Lems_derived_1 = NAN;
		// common read-only namespace? HHRate BaseTau 1 for Fixed channel 3
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_derived_0;
		float *Lems_assigned_4 = &Lems_derived_1;
		// compute derived HHRate BaseTau 1 for Fixed channel 3
		Lems_derived_0 = ( ( *Lems_assigned_0/* voltage */ / *Lems_assigned_2/* voltage */ ) )/* unitless */;
		Lems_derived_1 = ( ( *Lems_assigned_1/* time */ * ( ( ( ( ( ( 20/* unitless */ * ( expf( ( ( ( ( *Lems_assigned_3/* unitless */ + 160/* unitless */ ) )/* unitless */ / 30/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ / ( ( 1/* unitless */ + ( expf( ( ( ( ( *Lems_assigned_3/* unitless */ + 84/* unitless */ ) )/* unitless */ / 7.2999999999999998/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ + 35/* unitless */ ) )/* unitless */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 1 for Fixed channel 3
		float Lems_exposure_t = Lems_derived_1;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[38] / (1 + exp( (local_constants[39] - Vcomp ) / local_constants[40] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[4] = inf;
		}else{
			local_stateNext[4] = local_state[4] + dt * ( ( inf - local_state[4] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_1_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[41]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[42];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling ;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[43]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	// Current from synapses
	float I_synapses_total = 0;
	// Current from inputs
	float I_input_total = 0;
	I_internal = I_channels_total + I_input_total + I_synapses_total;
	if(initial_state){
		// initialize
		V_next[comp] = V[comp];
	}else{
		V_next[comp] = V[comp] + ( dt * ( I_internal ) / C[comp] );
	}
	}
	}
	// Internal Code for compartment type 1 end
	// Internal Code for compartment type 2
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[9];
	const long long Type_Compartments    = cell_const_table_i64_sizes [9];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	float Acomp = Area[comp];
	float Vcomp = V[comp];
	float I_internal = 0;
	// Ion flux sources
		float I_ion_0 = 0; //total ion current
		float Conc_ion_0_intra = 0; //ion concentration intra
		float Conc_ion_0_extra = 0; //ion concentration extra
	// Ion concentrations
	{
		float Iion = I_ion_0;
	float InitConcIntra = local_constants[0];
	float InitConcExtra = local_constants[1];
	Conc_ion_0_intra = local_state[0];
	Conc_ion_0_extra = local_state[1];
	}
	const float Ca_concentration = Conc_ion_0_intra;
	const float Ca_concentration_extra = Conc_ion_0_extra;
	const float Ca2_concentration = 0;
	const float Ca2_concentration_extra = 0;
	// Current from ion channels
	float I_channels_total = 0;
	{
	float Vshift = 0;
		float Erev  = local_constants[5];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[2]; 
	// dynamics for channel 0 gate 0 
	{
		float q10 = local_constants[6];
		float alpha;
		{
		float r; // define exposure
		r = local_constants[7] / (1 + exp( (local_constants[8] - Vcomp ) / local_constants[9] ) );
		alpha = r;
		}
		float beta;
		{
		float r; // define exposure
		r = local_constants[10] * ( ( Vcomp == local_constants[11]) ? 1 : ( ( (Vcomp - local_constants[11] ) / local_constants[12] )  / (1 - exp( - (Vcomp - local_constants[11] ) / local_constants[12] ) ) ) );
		beta = r;
		}
		float tau;
		tau = 1 / ( alpha + beta );
		float inf;
		inf = alpha / ( alpha + beta );
		if(initial_state){
			local_stateNext[2] = inf;
		}else{
			local_stateNext[2] = local_state[2] + dt * ( ( inf - local_state[2] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[13]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;
		I_ion_0 += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[14];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[3]; 
	// dynamics for channel 1 gate 0 
	{
		float q10 = local_constants[15];
		float alpha;
		{
		float r; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		float Lems_requirement_1 = Ca_concentration;
		// fixed properties HHRate BaseRate 0 for Fixed channel 1
		float Lems_property_0 = local_constants[16];
		float Lems_property_1 = local_constants[17];
		// state variables HHRate BaseRate 0 for Fixed channel 1
		// declare derived variables HHRate BaseRate 0 for Fixed channel 1
		float Lems_derived_0 = NAN;
		float Lems_derived_1 = NAN;
		// common read-only namespace? HHRate BaseRate 0 for Fixed channel 1
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_requirement_1;
		float *Lems_assigned_2 = &Lems_property_0;
		float *Lems_assigned_3 = &Lems_property_1;
		float *Lems_assigned_4 = &Lems_derived_0;
		float *Lems_assigned_5 = &Lems_derived_1;
		// compute derived HHRate BaseRate 0 for Fixed channel 1
		Lems_derived_0 = ( ( ( ( 2.0000000000000002e-05/* unitless */ * *Lems_assigned_1/* concentration */ ) )/* concentration */ / *Lems_assigned_3/* concentration */ ) )/* unitless */;
		Lems_derived_1 = 0;		if( 0 );
		else if( ( ( *Lems_assigned_4/* unitless */ > 0.01/* unitless */ ) )/* unitless */ ){
			Lems_derived_1 = ( ( 0.01/* unitless */ / *Lems_assigned_2/* time */ ) )/* per_time */;
		}
		else if( ( ( *Lems_assigned_4/* unitless */ <= 0.01/* unitless */ ) )/* unitless */ ){
			Lems_derived_1 = ( ( *Lems_assigned_4/* unitless */ / *Lems_assigned_2/* time */ ) )/* per_time */;
		}
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseRate 0 for Fixed channel 1
		float Lems_exposure_r = Lems_derived_1;
		r = Lems_exposure_r;
		alpha = r;
		}
		float beta;
		{
		float r; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseRate 0 for Fixed channel 1
		float Lems_property_0 = local_constants[18];
		// state variables HHRate BaseRate 0 for Fixed channel 1
		// declare derived variables HHRate BaseRate 0 for Fixed channel 1
		float Lems_derived_0 = NAN;
		// common read-only namespace? HHRate BaseRate 0 for Fixed channel 1
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_derived_0;
		// compute derived HHRate BaseRate 0 for Fixed channel 1
		Lems_derived_0 = ( ( 0.014999999999999999/* unitless */ / *Lems_assigned_1/* time */ ) )/* per_time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseRate 0 for Fixed channel 1
		float Lems_exposure_r = Lems_derived_0;
		r = Lems_exposure_r;
		beta = r;
		}
		float tau;
		tau = 1 / ( alpha + beta );
		float inf;
		inf = alpha / ( alpha + beta );
		if(initial_state){
			local_stateNext[3] = inf;
		}else{
			local_stateNext[3] = local_state[3] + dt * ( ( inf - local_state[3] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[19]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[20];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[4]; 
	// dynamics for channel 2 gate 0 
	{
		float q10 = local_constants[21];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 0 for Fixed channel 2
		float Lems_property_0 = local_constants[22];
		float Lems_property_1 = local_constants[23];
		// state variables HHRate BaseTau 0 for Fixed channel 2
		// declare derived variables HHRate BaseTau 0 for Fixed channel 2
		float Lems_derived_0 = NAN;
		float Lems_derived_1 = NAN;
		// common read-only namespace? HHRate BaseTau 0 for Fixed channel 2
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_derived_0;
		float *Lems_assigned_4 = &Lems_derived_1;
		// compute derived HHRate BaseTau 0 for Fixed channel 2
		Lems_derived_0 = ( ( *Lems_assigned_0/* voltage */ / *Lems_assigned_2/* voltage */ ) )/* unitless */;
		Lems_derived_1 = ( ( *Lems_assigned_1/* time */ / ( ( ( expf( ( ( ( ( ( ( - 0.085999999999999993/* unitless */ ) )/* unitless */ * *Lems_assigned_3/* unitless */ ) )/* unitless */ - 14.6/* unitless */ ) )/* unitless */ ) )/* unitless */ + ( expf( ( ( ( ( 0.070000000000000007/* unitless */ * *Lems_assigned_3/* unitless */ ) )/* unitless */ - 1.8700000000000001/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 0 for Fixed channel 2
		float Lems_exposure_t = Lems_derived_1;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[24] / (1 + exp( (local_constants[25] - Vcomp ) / local_constants[26] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[4] = inf;
		}else{
			local_stateNext[4] = local_state[4] + dt * ( ( inf - local_state[4] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[27]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[28];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
		{
		float x; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		float Lems_requirement_1 = Ca_concentration;
		// fixed properties HHRate BaseInf 0 for Fixed channel 3
		float Lems_property_0 = local_constants[29];
		float Lems_property_1 = local_constants[30];
		// state variables HHRate BaseInf 0 for Fixed channel 3
		// declare derived variables HHRate BaseInf 0 for Fixed channel 3
		float Lems_derived_0 = NAN;
		float Lems_derived_1 = NAN;
		float Lems_derived_2 = NAN;
		// common read-only namespace? HHRate BaseInf 0 for Fixed channel 3
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_requirement_1;
		float *Lems_assigned_2 = &Lems_property_0;
		float *Lems_assigned_3 = &Lems_property_1;
		float *Lems_assigned_4 = &Lems_derived_0;
		float *Lems_assigned_5 = &Lems_derived_1;
		float *Lems_assigned_6 = &Lems_derived_2;
		// compute derived HHRate BaseInf 0 for Fixed channel 3
		Lems_derived_1 = ( ( *Lems_assigned_1/* concentration */ / *Lems_assigned_3/* concentration */ ) )/* unitless */;
		Lems_derived_2 = ( ( 1/* unitless */ / ( ( 1/* unitless */ + ( expf( ( ( ( ( 0.00036999999999999999/* unitless */ - *Lems_assigned_5/* unitless */ ) )/* unitless */ / 0.089999999999999997/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */ ) )/* unitless */;
		Lems_derived_0 = ( ( *Lems_assigned_0/* voltage */ / *Lems_assigned_2/* voltage */ ) )/* unitless */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseInf 0 for Fixed channel 3
		float Lems_exposure_x = Lems_derived_2;
		x = Lems_exposure_x;
		chan_gate_0_q = x;
		}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[31]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[32];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling ;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[33]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	// Current from synapses
	float I_synapses_total = 0;
	// Current from inputs
	float I_input_total = 0;
 float iCa = I_ion_0; //total ion current
	// Dynamics for ion 0 
	{
		float Iion = I_ion_0;
	float InitConcIntra = local_constants[0];
	float InitConcExtra = local_constants[1];
	float ion_charge = 2;
	float influx_rate = NAN;
	influx_rate = ( (iCa / Acomp) * local_constants[4] );
	if(initial_state){
		// initialize
			local_stateNext[0] = local_state[0];
			local_stateNext[1] = local_state[1];
	}else{
			float leak_rate = ( ( local_state[0] - local_constants[2] ) / local_constants[3] );
			local_stateNext[0] = local_state[0] + ( dt * ( influx_rate - leak_rate ) );
			if( local_stateNext[0] < 0 ) local_stateNext[0] = 0;
			local_stateNext[1] = local_state[1];
	}
	}
	I_internal = I_channels_total + I_input_total + I_synapses_total;
	if(initial_state){
		// initialize
		V_next[comp] = V[comp];
	}else{
		V_next[comp] = V[comp] + ( dt * ( I_internal ) / C[comp] );
	}
	}
	}
	// Internal Code for compartment type 2 end
	// Internal Code for compartment type 3
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[10];
	const long long Type_Compartments    = cell_const_table_i64_sizes [10];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	float Acomp = Area[comp];
	float Vcomp = V[comp];
	float I_internal = 0;
	// Ion flux sources
	// Ion concentrations
	const float Ca_concentration = 0;
	const float Ca_concentration_extra = 0;
	const float Ca2_concentration = 0;
	const float Ca2_concentration_extra = 0;
	// Current from ion channels
	float I_channels_total = 0;
	{
	float Vshift = 0;
		float Erev  = local_constants[0];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
		{
		float x; // define exposure
		x = local_constants[1] / (1 + exp( (local_constants[2] - Vcomp ) / local_constants[3] ) );
		chan_gate_0_q = x;
		}
	float chan_gate_1_q; 
	chan_gate_1_q = local_state[0]; 
	// dynamics for channel 0 gate 1 
	{
		float q10 = local_constants[4];
		float tau;
		{
		float t; // define exposure
		// LEMS component
		float Lems_requirement_0 = Vcomp;
		// fixed properties HHRate BaseTau 1 for Fixed channel 0
		float Lems_property_0 = local_constants[5];
		float Lems_property_1 = local_constants[6];
		float Lems_property_2 = local_constants[7];
		// state variables HHRate BaseTau 1 for Fixed channel 0
		// declare derived variables HHRate BaseTau 1 for Fixed channel 0
		float Lems_derived_0 = NAN;
		// common read-only namespace? HHRate BaseTau 1 for Fixed channel 0
		float *Lems_assigned_0 = &Lems_requirement_0;
		float *Lems_assigned_1 = &Lems_property_0;
		float *Lems_assigned_2 = &Lems_property_1;
		float *Lems_assigned_3 = &Lems_property_2;
		float *Lems_assigned_4 = &Lems_derived_0;
		// compute derived HHRate BaseTau 1 for Fixed channel 0
		Lems_derived_0 = ( ( *Lems_assigned_1/* time */ * ( expf( ( ( ( ( *Lems_assigned_0/* voltage */ - *Lems_assigned_2/* voltage */ ) )/* voltage */ / *Lems_assigned_3/* voltage */ ) )/* unitless */ ) )/* unitless */ ) )/* time */;
		// integrate inline
		if(initial_state){
			// initialization
		}else{
			// dynamics
			// (highest up is lowest priority)
			// time derivatives
		// conditional updates, during simulation
		}
		// expose inline
		// exposures HHRate BaseTau 1 for Fixed channel 0
		float Lems_exposure_t = Lems_derived_0;
		t = Lems_exposure_t;
		tau = t;
		}
		float inf;
		{
		float x; // define exposure
		x = local_constants[8] / (1 + exp( (local_constants[9] - Vcomp ) / local_constants[10] ) );
		inf = x;
		}
		if(initial_state){
			local_stateNext[0] = inf;
		}else{
			local_stateNext[0] = local_state[0] + dt * ( ( inf - local_state[0] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_1_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[11]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[12];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float chan_gate_0_q; 
	chan_gate_0_q = local_state[1]; 
	// dynamics for channel 1 gate 0 
	{
		float q10 = local_constants[13];
		float alpha;
		{
		float r; // define exposure
		r = local_constants[14] * ( ( Vcomp == local_constants[15]) ? 1 : ( ( (Vcomp - local_constants[15] ) / local_constants[16] )  / (1 - exp( - (Vcomp - local_constants[15] ) / local_constants[16] ) ) ) );
		alpha = r;
		}
		float beta;
		{
		float r; // define exposure
		r = local_constants[17] * exp( (Vcomp - local_constants[18] ) / local_constants[19] );
		beta = r;
		}
		float tau;
		tau = 1 / ( alpha + beta );
		float inf;
		inf = alpha / ( alpha + beta );
		if(initial_state){
			local_stateNext[1] = inf;
		}else{
			local_stateNext[1] = local_state[1] + dt * ( ( inf - local_state[1] ) / tau ) * q10 ;
		}
	}
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling  * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q * chan_gate_0_q;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[20]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	{
	float Vshift = 0;
		float Erev  = local_constants[21];
	float ChannelOpenFraction = NAN;
	float ChannelConductance = NAN;
		float rateScale = 1;
	float conductance_scaling = 1;
		ChannelOpenFraction = conductance_scaling ;
	float I_chan = NAN;
	float iDensity = NAN;
		float Gbase = local_constants[22]; // conductivity
		float Gscaled  = Gbase * ChannelOpenFraction;
		iDensity = Gscaled * (Erev - Vcomp) * 1e-5f;
	I_chan = iDensity * Acomp;
		I_channels_total += I_chan;

	}

	// Current from synapses
	float I_synapses_total = 0;
	// Current from inputs
	float I_input_total = 0;
	I_internal = I_channels_total + I_input_total + I_synapses_total;
	if(initial_state){
		// initialize
		V_next[comp] = V[comp];
	}else{
		V_next[comp] = V[comp] + ( dt * ( I_internal ) / C[comp] );
	}
	}
	}
	// Internal Code for compartment type 3 end
	{
	const long long Compartments = cell_state_table_f32_sizes[0]; //same for all parallel arrays
	const Table_I64 Order  = cell_const_table_i64_arrays[11];
	const Table_I64 Parent = cell_const_table_i64_arrays[12];
	const Table_F32 DperT  = cell_const_table_f32_arrays[4];
	Table_F32 D = cell_state_table_f32_arrays[0];
	for(long long comp_seq = 0; comp_seq < Compartments; comp_seq++){
		D[comp_seq] = 1 + DperT[comp_seq] * dt ;
	}
	for( long long comp_seq = 0; comp_seq < Compartments - 1; comp_seq++ ){
		long long i = Order[comp_seq];
		long long j = Parent[i];
		long long idx = ( ( i > j ) ? i : j );
		float R = R_Axial[idx];
		float Ui = - dt/( R * C[i]) ;
		float Uj = - dt/( R * C[j]) ;
		float Li = Uj;
		float ratio = Li/D[i];
		D[j] -= ratio * Ui;
		V_next[j] -= ratio * V_next[i];
	}
	long long i = Order[ Compartments - 1 ];
	V_next[i] = V_next[i] / D[i];
	for( long long comp_seq = Compartments - 2; comp_seq >= 0 ; comp_seq-- ){
		long long i = Order[comp_seq];
		long long j = Parent[i];
		long long idx = ( ( i > j ) ? i : j );
		float R = R_Axial[idx];
		float Ui = - dt/( R * C[i]) ;
		V_next[i] = ( V_next[i] - Ui * V_next[j] ) / D[i];
	}
	}
	// PostUpdate Code for compartment type 0
	// Internal Code for compartment type 0
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[7];
	const long long Type_Compartments    = cell_const_table_i64_sizes [7];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	}
	}
	// Internal Code for compartment type 0 end
	// PostUpdate Code for compartment type 1
	// Internal Code for compartment type 1
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[8];
	const long long Type_Compartments    = cell_const_table_i64_sizes [8];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	}
	}
	// Internal Code for compartment type 1 end
	// PostUpdate Code for compartment type 2
	// Internal Code for compartment type 2
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[9];
	const long long Type_Compartments    = cell_const_table_i64_sizes [9];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	}
	}
	// Internal Code for compartment type 2 end
	// PostUpdate Code for compartment type 3
	// Internal Code for compartment type 3
	{
	const Table_I64 Comp_List    = cell_const_table_i64_arrays[10];
	const long long Type_Compartments    = cell_const_table_i64_sizes [10];
	for( long long CompIdx = 0; CompIdx < Type_Compartments; CompIdx++ ){
		int comp = (int) Comp_List[CompIdx];
		const long long const_comp_index      = Comp_Coff   [comp];
		const long long state_comp_index      = Comp_Soff   [comp];
		const long long table_cf32_comp_index = Comp_CF32off[comp];
		const long long table_ci64_comp_index = Comp_CI64off[comp];
		const long long table_sf32_comp_index = Comp_SF32off[comp];
		const long long table_si64_comp_index = Comp_SI64off[comp];
		const long long rng_offset            = Comp_Roff   [comp];
		
	const float *comp_constants = cell_constants + const_comp_index;
	const float *comp_state     = cell_state     + state_comp_index;
	      float *comp_stateNext = cell_stateNext + state_comp_index;
	
		const long long *comp_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_comp_index;
		const Table_F32 *comp_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_comp_index;
		const long long *comp_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_comp_index;
		const Table_I64 *comp_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_comp_index;
		const long long *comp_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_comp_index;
		const Table_F32 *comp_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_comp_index;
		      Table_F32 *comp_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_comp_index;
		const long long *comp_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_comp_index;
		      Table_I64 *comp_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_comp_index;
		      Table_I64 *comp_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_comp_index;
	const long long const_local_index = const_comp_index;
	const long long state_local_index = state_comp_index;
	const long long table_cf32_local_index = table_cf32_comp_index;
	const long long table_ci64_local_index = table_ci64_comp_index;
	const long long table_sf32_local_index = table_sf32_comp_index;
	const long long table_si64_local_index = table_si64_comp_index;
	
	const float *local_constants = cell_constants + const_local_index;
	const float *local_state     = cell_state     + state_local_index;
	      float *local_stateNext = cell_stateNext + state_local_index;
	
		const long long *local_const_table_f32_sizes      = cell_const_table_f32_sizes      + table_cf32_local_index;
		const Table_F32 *local_const_table_f32_arrays     = cell_const_table_f32_arrays     + table_cf32_local_index;
		const long long *local_const_table_i64_sizes      = cell_const_table_i64_sizes      + table_ci64_local_index;
		const Table_I64 *local_const_table_i64_arrays     = cell_const_table_i64_arrays     + table_ci64_local_index;
		const long long *local_state_table_f32_sizes      = cell_state_table_f32_sizes      + table_sf32_local_index;
		const Table_F32 *local_state_table_f32_arrays     = cell_state_table_f32_arrays     + table_sf32_local_index;
		      Table_F32 *local_stateNext_table_f32_arrays = cell_stateNext_table_f32_arrays + table_sf32_local_index;
		const long long *local_state_table_i64_sizes      = cell_state_table_i64_sizes      + table_si64_local_index;
		      Table_I64 *local_state_table_i64_arrays     = cell_state_table_i64_arrays     + table_si64_local_index;
		      Table_I64 *local_stateNext_table_i64_arrays = cell_stateNext_table_i64_arrays + table_si64_local_index;
	}
	}
	// Internal Code for compartment type 3 end
}
// Generated code block END
