<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LFRHU diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRHU_LFRHU_SEM_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">Couette_FRHU_LFRHU_SEM_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-09">0.000561628</value>
            <value variable="rhou" tolerance="1e-04">68.0654</value>
            <value variable="rhov" tolerance="1e-06">0.206241</value>
            <value variable="rhow" tolerance="1e-10">1.32025e-05</value>
            <value variable="E" tolerance="1e-01">24776.9</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-09">0.00139967</value>
            <value variable="rhou" tolerance="1e-04">83.3517</value>
            <value variable="rhov" tolerance="1e-06">0.50519</value>
            <value variable="rhow" tolerance="1e-10">3.11554e-05</value>
            <value variable="E" tolerance="1e-01">18953</value>
        </metric>
    </metrics>
</test>


