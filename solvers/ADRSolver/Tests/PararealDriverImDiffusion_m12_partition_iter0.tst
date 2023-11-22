<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion, order 1, P=12 </description>
    <executable>ADRSolver</executable>
    <parameters>--npt 5 PararealDriverImDiffusion_m12_iter0_xml</parameters>
    <processes>10</processes>
    <files>
        <file description="Session File"> PararealDriverImDiffusion_m12_iter0_xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">0.00887119</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">0.0117097</value>
        </metric>
    </metrics>
</test>
