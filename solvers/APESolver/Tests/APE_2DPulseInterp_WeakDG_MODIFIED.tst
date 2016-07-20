<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>APESolver</executable>
    <parameters>APE_2DPulseInterp_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DPulseInterp_WeakDG_MODIFIED.xml</file>
        <file description="Pts File">APE_2DPulseInterp_baseflow_0.00000000E+00.pts</file>
        <file description="Pts File">APE_2DPulseInterp_baseflow_1.00000000E-05.pts</file>
        <file description="Pts File">APE_2DPulseInterp_baseflow_2.00000000E-05.pts</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-6">22.1182</value>
            <value variable="u" tolerance="1e-6">0.00207907</value>
            <value variable="v" tolerance="1e-6">0.00207907</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="5e-4">99.6754</value>
            <value variable="u" tolerance="1e-6">0.0079465</value>
            <value variable="v" tolerance="1e-6">0.0079465</value>
        </metric>
    </metrics>
</test>
