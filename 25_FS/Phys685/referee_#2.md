# [The Data Acquisition System of the LZ Dark Matter Detector: FADR](https://arxiv.org/pdf/2405.14732)
# Andrew Koren

Overview (one-line outline): The paper describes how the 745 PMTs used in the Xenon TPC are controlled, in particular the hardware control as well as the methods used to maximized data collection on events while minimze compute and memory waste on noise.

## Outline

Introduction:
-   LUX-ZEPLIN (LZ) is the follow up to the LUX (Large Underground Xenon) detector, featuring a liquid xenon tank inside of a water tank, with a titanium cryostat and a skin of gadolinium liquid scintillator between them.
-   PMTs are arranged above and below the Xenon tank (TPC), with skin PMTs facing the scintillator (Skin) and outer PMTs facing the water tank (OD). 
-   All detection components are controlled down to the firmware level by a single unified data acquisition system, FADR.

Analog-to-Digital Conversion & Hardware:
-   PMT power and output signal travel through custom vacuum flanges, with amplifiers in the breakout boxes on the outside.
-   Amplified signals from TPC, Skin and OD are routed to FADR for digitization
-   FADR uses 47 32-channel digitizers, one of which is dedicated to complimentary sensors, and one dedicated to  DAQScope
-   Each digitizer has on-board waveform memory, RAM, Arm processor, and gigabit ethernet
-   Alongside digital output, each board outputs diagnostics and clock time can be set by outside signal for synchronization
-   Each board outputs digital signal to a Data Sparsifier/Master which determines if an event of interest has occured and report to the Data Aquisition System (DAQ) Master. From there, the DAQ master determines the timeframe of the event occurence and orders all relevant data from the digitizes be stored for offline analysis.

FADR Design & Firmware:
-   The S1/S2 signal time difference is up to 700µs, so the system is designed to buffer this delay while limiting calibration rates to 150 Hz to prevent event pile-up.
-   To reduce data volume by a factor of 50, the firmware utilizes Pulse Only Digitization (POD), which only records waveform segments that cross a dynamic threshold.
-   When ordered by the DAQ Master, Data Extractors pull waveform segments from the digitizers and transmit them to Data Collectors over Ethernet via User datagram Protocol (UDP)
-   Event Builders compile data from the Collectors into the final product, complete event files, which are transferred from local 72TB disk space to RAID arrays at Sanford before transfer to NERSC data center at LBNL.
-   Time-stamping is run by a global 100-MHz clock distributed via HDMI cables, synchronized with a surface GPS antenna for absolute timestamping.
-   The firmware implements real-time Zero Suppression by calculating a rolling baseline and only storing samples when the signal deviates significantly. Multiple triggers are merged if their buffer time regions overlap. 
-   The last step carried out by digitizer is Digital Finite-Impulse-Response (FIR) filtering. This analyzes waveforms in real-time to identify S1 (narrow) and S2 (wide) pulses, sending this information to the Data Sparsifier to decide if an event should be kept. 

Event characterization:
-   During standard WIMP searches, the system uses a logical OR of three triggers: the TPC S2 trigger, the Skin S1 trigger, and the Outer Detector S1 trigger.
-   To ensure the detector is sensitive to S2-only events in the xenon, the TPC trigger logic is specifically tuned to the broader S2 pulse shape.
-   External triggers are employed for health monitoring, including a random trigger for unbiased noise sampling and a GPS trigger that runs at exactly 1 Hz to verify time.
-   Specific calibration triggers are used for maintenance, such as LED triggers for PMT response monitoring and Deuterium-Deuterium (DD) triggers for neutron generation.
-   The size of a standard WIMP search event is approximately 0.9 MB, most of which is due to the width of S2 signals
-   Calibration event sizes vary significantly, ranging from small 0.06 MB files for LED tests to large 2.1 MB files for neutron calibrations caused by multiple scattering interactions.

Extended Functionality:
-   Firmware allows for arbitrary waveform injection at the digitizer front-end, enabling the testing of trigger logic and data processing paths using synthetic pulses mixed with real noise.
-   DAQScope provides a virtual oscilloscope capability, allowing operators to monitor individual or summed channels remotely via video streams without disrupting physical connections.
-   A digital sum feature aggregates waveforms from multiple digitizers into a single data stream, which the Data Sparsifier Master uses to make global trigger decisions based on total pulse area.

Performance of all systems is also measured with a variety of criteria.

## Referee

### Introduction

The LUX-ZEPLIN (LZ) experiment, deployed at the Sanford Underground Research Facility (SURF), represents the next generation of dark matter detection, succeeding the Large Underground Xenon (LUX) detector. At its core is a dual-phase Time Projection Chamber (TPC) containing liquid xenon, housed within a titanium cryostat. To minimize background interference, this cryostat is surrounded by a skin detector (skin) composed of gadolinium-loaded liquid scintillator (GdLS), and the entire assembly is submerged inside a large water tank acting as an Outer Detector (OD).

The TPC is monitored by 494 3-inch PMTs arranged in arrays above and below the xenon volume. The Skin detector utilizes 93 1-inch and 20 2-inch PMTs looking downwards, along with 18 2-inch PMTs in the dome region. The Outer Detector is monitored by 120 8-inch PMTs mounted on stainless steel ladders within the water tank. All detection components, from the TPC to the OD, are readout and controlled down to the firmware level by a single unified system: the FPGA-based Architecture for Data acquisition and Realtime monitoring (FADR).

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_LZDetector.png"
         alt="Figure 1: Schematic of the LZ detector with its main components indicated, including the TPC, cryostat, and water tank."
         width="750">
    <figcaption>Figure 1: Schematic of LZ detector. Below are a list of components.
    <li> TPC (1) and water tank (2) </li> 
    <li> top PMTs (3) and bottom PMTs (4) internal to TPC </li>
    <li> titanium cryostat (5)  </li> 
    <li> liquid scintillator (6) and skin PMTs (7)(8) </li>
    <li> PMT signals are routed through the top (9) and bottom (10) conduits </li>
    </figcaption>
</figure>

The system is designed to detect two distinct types of signals generated by particle interactions: the prompt scintillation signal (S1) and the ionization signal (S2). The S1 signal is extremely fast, with a full width at half maximum (FWHM) of less than 100 ns. The S2 signal, caused by electroluminescence as electrons are extracted into the gas phase, occurs up to 700 µs later—depending on the depth of the interaction—and is significantly broader, spanning several microseconds. This dual-signal characteristic allows for precise 3D position reconstruction: the S2 light distribution on the top array indicates horizontal position, while the time delay between S1 and S2 reveals the vertical depth.

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_Xe_TPC.png"
         alt="Figure 2: Principle of operation of a dual-phase xenon detector, showing S1 and S2 signal generation and PMT arrays."
         width="750">
    <figcaption>Figure 2: Illustration of event in LZ. L2 electrons drift in the TPC before being observed via electroluminescence. The top PMTs are colored to match an L2 event's amplitudes.</figcaption>
</figure>

### Analog-to-Digital Conversion & Hardware

PMT voltage signals originating from within the cryostat are routed via thin coaxial cables (Axon) to custom-made vacuum flanges. Each flange is equipped with four DB-25 feedthroughs to allow data transfer between the pressure/temperature gap.

Analog signals are immediately processed by custom 8-channel dual-gain amplifiers mounted directly onto the signal flanges in 5-card mini crates. These amplifiers split the incoming PMT signals into a slow (60ns) high-gain channel and a fast (30ns) low-gain channel. The amplified signals are then transmitted to the FADR electronics racks via LMR-100A cables for the xenon space and RG316/U cables for the OD.

Powering these PMTs requires distinct high-voltage (HV) strategies. The TPC and Skin PMTs operate at a negative high voltage, supplied by WIENER MPOD (awsome product name) modules through filters mounted directly on the breakout box flanges. Outer Detector PMTs operate at a positive high voltage, with HV filters installed directly to the power supplies, rather than on breakout boxes.

The analog signals are sent to 47 32-channel digitizers (DDC-32s) and 24 logic boards installed in 6U VME crates. Of these digitizers, 45 are dedicated to processing the 1,359 PMT channels, one monitors fast sensors for environmental diagnostics, and one is dedicated to the DAQScope monitoring system. The DDC-32 is a custom board developed by Skutek, capable of sampling at 100 MHz with 14-bit resolution over a 2-V range (-1.8 V to +0.2 V).

The processing power of the DDC-32 is provided by an onboard Xilinx Kintex-7 Field Programmable Gate Array (FPGA). For control / diagnostics, each digitizer uses an Arm Cortex-8 (AM3358) and 512 MB of DDR3 RAM to run a Linux operating system, allowing for flexible operation and network communication.

Connectivity is handled through multiple interfaces. Waveform data and trigger information are transmitted via HDMI cables using Low-Voltage Differential Signaling (LVDS) to logic boards. Additionally, each board features Gigabit Ethernet for data readout and RS-232 for local debugging. The logic boards aggregate data from the digitizers to the Data Sparsifier Master, which makes the final determination on whether an event of interest has occurred. The Master then communicates with the central DAQ Master to orchestrate the readout of relevant time windows from all digitizers for offline analysis.

### FADR Design & Firmware
The FADR firmware is engineered to handle the specific timing challenges of the TPC, specifically the 700 µs drift time between S1 and S2 signals. To manage data volume efficiently, the system utilizes Pulse Only Digitization (POD). Rather than recording continuous empty baselines, the system only stores waveform segments that cross a threshold. This reduces data volume by a factor of approximately 50, even before the Data Sparsifier Master throws out false hits. The POD logic includes 32 samples of "pre-trigger" and 32 samples of "post-trigger" data around the pulse to ensure the baseline and pulse tails are fully captured.

The data flow is a multi-stage process. When the DAQ Master triggers an event readout, Data Extractors pull the relevant PODs from the DDC-32 circular buffers via LVDS links. These extractors then transmit the data to Data Collectors using User Datagram Protocol (UDP) over Gigabit Ethernet, a different scheme than the TCP used by most computer's ethernet connection. The Data Collectors temporarily store the data on local solid-state drives before Event Builders compile the fragments into complete event files. These files are moved to a 72 TB surface RAID array and subsequently transferred to the NERSC data center at LBNL for permanent storage.

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_DAQTopLevel.png"
         alt="Figure 7: Schematic view of the operation of FADR, showing circular buffers and digital filters on the FPGA."
         width="750">
    <figcaption>Figure 7: Schematic view of the operation of FADR. The Circular Buffer is made up of two buffers used to load data into/out of the digitizer without pause. I'm not sure what a Summer Selector is, maybe that is the "spy" system? Also, The Data Extractor is controlled by the DAQ DB, but there's no reference to that in this diagram. Otherwise, it sums up the past two sections pretty well.</figcaption>
</figure>

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_FADR_Architecture.png"
         alt="Figure 8: Diagram of the architecture of FADR, showing groups of DDC-32s, Data Extractors, and the DAQ Master."
         width="750">
    <figcaption>Figure 8: A more zoom-in to the event-builder portion of the previous diagram. The top-right DDC-32 is for non-PMT sensors, such as accoustic data, which is bundled into the event report.</figcaption>
</figure>

System timing is maintained by a global 100-MHz clock distributed to all boards via the differential pairs in the HDMI cables. For absolute timestamping, the system is synchronized with a GPS receiver located on the surface, which provides a 1-pulse-per-second (PPS) signal. This ensures that every sample recorded by the FADR system can be correlated with Coordinated Universal Time (UTC) to a precision of 15 ns.

A critical component of the firmware is the real-time event selection logic. The FPGA applies digital Finite-Impulse-Response (FIR) filters to the incoming data streams to distinguish between S1 and S2 pulses. The S1 filter looks for narrow, prompt signals, while the S2 filter is tuned for wider pulses. This multiplicity and pulse-area information is sent to a Data Sparsifier Master, which aggregates inputs from all digitizers to determine if a global trigger condition (logical OR of TPC S2, Skin S1, or OD S1) has been met.

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_S1S2Filters.png"
         alt="Figure 14: The FIR filters used to process the incoming waveforms to look for S1 and S2 signals."
         width="750">
    <figcaption>Figure 14: The FIR filters used to process the incoming waveforms to look for S1 and S2 signals. The diagram shows how side lobes are weighted negatively compared to the main lobes, with S1 matching the tighter signal (N=160ns integration time) and S2 the wider, deeper signal (M=5.12µs integration time). This separates the two and helps prevent slow signal fluctuations from breaking the dynamic threshold.
    </figcaption>
</figure>

### Event Characterization and Performance
During standard WIMP search operations, the event selection logic is a logical OR of three primary triggers: the TPC S2 trigger, the Skin S1 trigger, and the Outer Detector S1 trigger. The TPC trigger is unique in that it is specifically tuned to the broader shape of the electroluminescence signal, ensuring the detector remains sensitive to S2-only events which may lack a detectable S1 pulse. The Skin and OD triggers focus on the prompt S1-like signals characteristic of gamma or neutron scattering in those volumes.

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_eventDefinition.png"
         alt="Figure 9: Event definition showing the pre-event and post-event windows."
         width="750">
    <figcaption>Figure 9: Event definition showing the pre-event and post-event windows. The DAQ Master determines this interval algorithmically, and orders the Data Extractor logic boards to collect and send high-resolution data from the digitizers to the Data Collectors, alongside various hardware reports / basic signal analysis. </figcaption>
</figure>


Beyond physics triggers, the system employs external triggers for health and consistency monitoring. A random trigger operates at approximately 4 Hz to provide an unbiased sampling of the noise environment. Simultaneously, a GPS trigger forces an acquisition exactly once per second (1 Hz). This GPS trigger is crucial for verifying the system's "livetime"—the actual duration the detector is active and capable of recording data.

<!-- Calibration triggers are essential for detector maintenance. LED triggers are used to flash light into the detector to monitor PMT gains and response, producing very small event files (~0.06 MB). Conversely, Deuterium-Deuterium (DD) neutron generator triggers create complex scattering events throughout the xenon and veto regions, resulting in much larger event files (~2.1 MB). For comparison, a standard WIMP search event has a file size of approximately 0.9 MB, a volume dominated by the long digitization window required for the wide S2 signals. -->

As mentioned sporatically in this review, FADR includes features that aid in debugging and monitoring. The firmware supports Arbitrary Waveform Injection at the front-end of the digitizers. This allows operators to inject synthetic pulses into the signal chain after the ADCs, allowing for real-world testing of triggers and logic boards rather than relying on simulated devices.

Another key feature is DAQScope, which functions as a virtual oscilloscope. It allows operators to select specific analog channels to be routed to a "spy" output. This output is converted back to analog and digitized by a video encoder, enabling remote viewing of real-time waveforms via a video stream without physically disturbing the cabling or interrupting data acquisition.

Also, The Data Sparsifier Master is able to make inclusion descisions from multiple digitzers using a Digital Sum capability. This feature aggregates the digitized waveforms from multiple DDC-32 boards into a single data stream. The Data Sparsifier Master then uses the total pulse area across the entire detector, rather than relying solely on individual channel thresholds, to determine an event has triggered.

The performance of FADR is continuously monitored against strict criteria. The intrinsic noise of the DDC-32 channels is measured to be between 2.2 and 2.4 ADC counts (RMS), ensuring a high signal-to-noise ratio for single photoelectron detection. The integral nonlinearity of the system is maintained below ±3 ADC counts, with differential nonlinearity much lower.

In terms of data throughput metrics, the system is highly capable. Data acquisition rates get up to 1200 MB/s (~0.9MB for a single event). Livetime (uptime where it can record the full S1/S2 gap) is a critical metric for a rare-event search. During WIMP search modes, the livetime is approximately 99.5%. When a 2 ms post-event holdoff is applied to prevent re-triggering on S2 tails, the effective livetime is reduced slightly to 95.9%, which is well within acceptable parameters for the experiment.

## Conclusion

Overall, I recommend accepting the paper due to its completeness in discussion of design and implementation. The paper does well to establish a background for those outside the exact technical specialization of data acquisition, while still explaining and revealing every technical detail that went into overcoming the S1/S2 timegap, hardware choices, and debugging tools that were probably developed to prevent headaches during the LUX experiment. The layout is straightforward.

Publish as is.

<!-- ----

### Introduction and Detector Overview

The LUX-ZEPLIN (LZ) experiment, deployed at the Sanford Underground Research Facility (SURF), represents the next generation of dark matter detection, succeeding the Large Underground Xenon (LUX) detector. At its core is a dual-phase Time Projection Chamber (TPC) containing liquid xenon, housed within a titanium cryostat. To minimize background interference, this cryostat is surrounded by a skin detector (skin) composed of gadolinium-loaded liquid scintillator (GdLS), and the entire assembly is submerged inside a large water tank acting as an Outer Detector (OD). The detection system utilizes arrays of photomultiplier tubes (PMTs) to capture signals: 253 3-inch PMTs are positioned in the gas phase above the xenon, and 241 are immersed in the liquid below. Additionally, the Skin detector utilizes 93 1-inch and 20 2-inch PMTs monitoring the scintillator layer, while the OD employs 120 8-inch PMTs to identify neutron-scattering events in the water tank.


The system is designed to detect two distinct types of signals generated by particle interactions: the prompt scintillation signal (S1) and the ionization signal (S2). The S1 signal is extremely fast, with a full width at half maximum (FWHM) of less than 100 ns. The S2 signal, caused by electroluminescence as electrons are extracted into the gas phase, occurs up to 700 µs later—depending on the depth of the interaction—and is significantly broader, spanning several microseconds. This dual-signal characteristic allows for precise 3D position reconstruction: the S2 light distribution on the top array indicates horizontal position, while the time delay between S1 and S2 reveals the vertical depth.

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_Xe_TPC.png"
         alt="Figure 2: Principle of operation of a dual-phase xenon detector, showing S1 and S2 signal generation and PMT arrays."
         width="750">
    <figcaption>Figure 2: Illustration of event in LZ. L2 electrons are shot directionally and then observed via electroluminescence.</figcaption>
</figure>

### Analog-to-Digital Conversion & Hardware

The analog signal chain begins with PMT signals exiting the cryostat through custom vacuum flanges equipped with DB-25 feedthroughs, located on breakout boxes on the water tank's exterior. To maximize dynamic range, the signals are processed by 8-channel dual-gain amplifiers mounted directly on these flanges. These amplifiers split the signal into high-gain (area gain of 40) and low-gain (area gain of 4) channels. The conditioned signals are then routed via LMR coaxial cables to the FADR system for digitization.

The input to the FADR system is a bank of 47 32-channel digitizers (DDC-32s). These boards are built around the Xilinx Kintex-7 Field Programmable Gate Array (FPGA) and digitize signals at 100 MHz with 14-bit resolution over a 2-V range (-1.8 V to +0.2 V). Each DDC-32 computes with 512 MB of DDR3 RAM and an Arm Cortex-8 (AM3358) processor running Linux for control and diagnostics. The system operates synchronously via a global 100-MHz clock distributed via HDMI cables, which is aligned with a 1-pulse-per-second GPS signal for absolute timestamping.

### FADR Architecture and Firmware

The FADR firmware is engineered to handle the specific timing challenges of the TPC, specifically the 700 µs drift time between S1 and S2 signals. To manage data volume efficiently, the system utilizes Pulse Only Digitization (POD). Rather than recording the continuous baseline, the FPGA calculates a rolling average and only stores waveform segments that cross a dynamic threshold. This reduces data volume by a factor of approximately 50, even before the Data Sparsifier Master throws out false hits. The firmware implements a dual-buffer system where samples are stored in circular buffers; when an event is triggered, Data Extractors retrieve specific waveform segments from the digitizers via LVDS links and transmit them to Data Collectors over Ethernet using the UDP protocol. This point-to-point UDP topology eliminates TCP overhead, allowing for a sustained total data acquisition rate of 1,200 MB/s.

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_FADR_Architecture.png"
         alt="Figure 7: Schematic view of the operation of FADR, showing circular buffers and digital filters on the FPGA."
         width="750">
    <figcaption>Figure 7: Schematic view of the operation of FADR. I'm not sure what a Summer Selector is.</figcaption>
</figure>


<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_FADR_Architecture.png"
         alt="Figure 8: Diagram of the architecture of FADR, showing groups of DDC-32s, Data Extractors, and the DAQ Master."
         width="750">
    <figcaption>Figure 8: Illustration of event in LZ. L1 direct photon emission happens near instantly, while electrons are shot directionally and then observed via electroluminescence.</figcaption>
</figure>

A critical component of the firmware is the real-time event selection logic. The FPGA applies digital Finite-Impulse-Response (FIR) filters to the incoming data streams to distinguish between S1 and S2 pulses. The S1 filter looks for narrow, prompt signals, while the S2 filter is tuned for wider pulses. This multiplicity and pulse-area information is sent to a Data Sparsifier Master, which aggregates inputs from all digitizers to determine if a global trigger condition (logical OR of TPC S2, Skin S1, or OD S1) has been met.

<figure>
    <img src="https://arxiv.org/html/2405.14732v3/extracted/5796747/Figs_S1S2Filters.png"
         alt="Figure 14: The FIR filters used to process the incoming waveforms to look for S1 and S2 signals."
         width="750">
    <figcaption>Figure 8: Illustration of event in LZ. L1 direct photon emission happens near instantly, while electrons are shot directionally and then observed via electroluminescence.</figcaption>
</figure>


### Event Characterization and Performance

Once an event is identified, the DAQ Master orders the extraction of data within a specific time window, defined by pre-event and post-event parameters. These windows are configurable; for example, during standard WIMP searches, the window captures the entire drift time to ensure both S1 and S2 signals are paired correctly.

The system supports various trigger modes to accommodate different operational needs. During calibration, "random" triggers provide unbiased noise sampling, while LED triggers (running at up to 4 kHz) monitor PMT health. For internal physics calibrations using neutron sources (like the Deuterium-Deuterium generator), the rate is throttled to 150 Hz to prevent pile-up, given the long electron drift time. The resulting data footprint varies by mode: a standard WIMP search event is approximately 0.9 MB, largely dominated by the S2 pulses, whereas neutron calibration events can reach 2.1 MB due to multiple scattering interactions. The system also features extended functionality such as "DAQScope," which allows operators to view real-time video streams of waveform data from any channel via the "spy" outputs on the DDC-32s, facilitating remote diagnostics without interrupting data taking. -->