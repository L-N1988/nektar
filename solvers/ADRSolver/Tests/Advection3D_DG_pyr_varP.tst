<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, pyramid elements, variable polynomial orders </description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_DG_pyr_varP.xml</parameters>
    <files>
        <file description="Session and Mesh File">Advection3D_DG_pyr_varP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">0.000700559</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">0.0141906</value>
        </metric>
    </metrics>
</test>