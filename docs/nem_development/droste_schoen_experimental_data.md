# Droste & Schoen (1988) Experimental Data Summary

## Tank Specifications (Identical for all 3 experiments)

- **Gross volume**: 4.85 m³
- **Tank diameter**: 1.25 m
- **Tank length**: 4.3 m
- **Cylindrical part length**: 3.6 m
- **Orientation**: Horizontal
- **Head type**: Korbbogen-type (DIN 28013) - Torispherical, similar to DIN standard
- **Material**: StE 36 (fine-grained steel, yield strength 360 N/mm²)
- **Vessel density**: 7700 kg/m³
- **Heat capacity**: 500 J/kg·K

### Wall Thickness
- **Test 1**: Cylindrical = 5.9 mm, Heads = 6.6 mm
- **Tests 2 & 3**: Cylindrical = 6.4 mm, Heads = 6.8 mm

## Tank Equipment

- **Safety valve**: 1" with set pressure of 15.6 bar
- **Liquid discharge valve**: 3/4"
- **Filling valve**: 1 1/4"
- **Gas discharge valve**: 3/4"

## Filling Conditions

- **Fill level**: 50% for all tests = 1.425 m³ liquid propane
- **Liquid level (as fraction of diameter)**: 0.5

## Fire Conditions

- **Type**: Surrounding fuel oil pool fire (NOT full engulfment)
- **Configuration**: Steel troughs 30 cm from tank, 60 cm wide
- **Simulates**: Pedestal-mounted tank

## Experimental Results

### Test 1 (October 1982)
- **Ambient temperature**: 10°C
- **Initial propane temperature**: 10°C (283 K)
- **Initial propane pressure**: 5.5 bar (550,000 Pa)
- **Time to PSV opening**: 5'40" = 340 seconds
- **PSV opening pressure**: 16.4 bar
- **Time to rupture**: 12'00" = 720 seconds
- **Rupture pressure**: 24.5 bar
- **Liquid temperature at rupture**: 72°C (345 K)
- **Pressure rise from ignition to rupture**: 19 bar

### Test 2 (November 1983)
- **Ambient temperature**: 2°C
- **Initial propane temperature**: 37°C (310 K)
- **Initial propane pressure**: 13.5 bar (1,350,000 Pa)
- **Time to PSV opening**: 1'40" = 100 seconds
- **PSV opening pressure**: 17.3 bar
- **Time to rupture**: 7'20" = 440 seconds
- **Rupture pressure**: 39 bar
- **Liquid temperature at rupture**: 84-87°C (357-360 K)
- **Pressure rise from ignition to rupture**: 25.5 bar

### Test 3 (December 1983)
- **Ambient temperature**: -3°C
- **Initial propane temperature**: 26°C (299 K)
- **Initial propane pressure**: 9.8 bar (980,000 Pa)
- **Time to PSV opening**: 2'30" = 150 seconds
- **PSV opening pressure**: 16.0 bar
- **Time to rupture**: 9'00" = 540 seconds
- **Rupture pressure**: 30.5 bar
- **Liquid temperature at rupture**: 77-78°C (350-351 K)
- **Pressure rise from ignition to rupture**: 20.7 bar

## Key Observations

1. **Liquid temperature rise rate**: ~6.7-7.4°C/min after 1-2 min delay for full fire development
2. **Pressure rise rate**:
   - Test 1: 1.8 bar/min
   - Test 2: 4.1 bar/min
   - Test 3: 2.6 bar/min
3. **Time to rupture correlation**: tR = 13 - 0.154·TB (where TB is initial propane temperature in °C)
4. **PSV discharge capacity**: 64 m³/min (air at normal conditions)
5. **PSV effectiveness**: Not sufficient to prevent pressure rise in fire scenario

## Parameters for HydDown Simulation

### Fixed Parameters (Same for all experiments)
- Vessel volume: 4.85 m³
- Diameter: 1.25 m
- Length: 4.3 m
- Orientation: horizontal
- Type: "DIN" (closest to Korbbogen/DIN 28013)
- Fill level: 0.5 (50%)
- Fluid: propane
- PSV set pressure: 15.6 bar = 1,560,000 Pa
- PSV discharge coefficient: 0.975
- Blowdown: 0.25 (typical)
- Back pressure: 101,325 Pa (atmospheric)

### Adjustable Parameters (for tuning)
- **PSV diameter**: Start with API 526 "D" = 0.0111 m (to be tuned)
- **Fire type**: "api_pool" (60 kW/m²) or "scandpower_pool" (100 kW/m²) (to be tuned)
- **Fire scaling**: Start with 1.0, may need < 1.0 for partial engulfment (to be tuned)

### Variable Parameters (Different for each test)
- Initial temperature: Test1=283K, Test2=310K, Test3=299K
- Initial pressure: Test1=550kPa, Test2=1350kPa, Test3=980kPa
- Wall thickness: Test1=0.0059m, Tests2&3=0.0064m
- End time: Test1=800s, Test2=600s, Test3=700s

## Rupture Material

- Material: CS_360LT (corresponds to StE 36 with 360 N/mm² yield strength)
- Fire type for rupture: Should match heat_transfer fire type
