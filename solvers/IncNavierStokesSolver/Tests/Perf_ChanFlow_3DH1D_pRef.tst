<?xml version="1.0" encoding="utf-8"?>
<test runs="15">
    <description>Laminar Channel Flow 3D homogeneous 1D, p-Refinenement tag</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Perf_ChanFlow_3DH1D_pRef.xml</parameters>
    <files>
        <file description="Session File">Perf_ChanFlow_3DH1D_pRef.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.00596321</value>
            <value variable="v" tolerance="1e-8">0.000227709</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-7">0.0100519</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0.015886</value>
            <value variable="v" tolerance="1e-8">0.00109468</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-7">0.0671074</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="8e-1" hostname="42.debian-bullseye-performance-build-and-test">13.758</value>
        </metric>
    </metrics>
</test>
