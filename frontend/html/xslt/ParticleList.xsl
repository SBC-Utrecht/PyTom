<?xml version="1.0" encoding="UTF-8"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

	<xsl:template match="/ParticleList">
				
				<table class="particleList">
					<xsl:for-each select="Particle">
						
						<tr class="particleList">
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:value-of select="@Filename"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							<xsl:element name="b">
								P<xsl:value-of select="position()"/>
							</xsl:element> 
						</xsl:element>
							
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:apply-templates select="Rotation"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							Rotation
						</xsl:element>			
							
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:apply-templates select="Shift"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							Shift
						</xsl:element>	
							
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:apply-templates select="WedgeInfo"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							Wedge
						</xsl:element>
							
						</tr>
					</xsl:for-each>
				</table>
				
		
	</xsl:template>
	
	<xsl:template match="Rotation">Z1 <xsl:value-of select="@Z1"/> X <xsl:value-of select="@X"/> Z2 <xsl:value-of select="@Z2"/></xsl:template>
	
	<xsl:template match="Shift">X <xsl:value-of select="@X"/> Y <xsl:value-of select="@Y"/> Z <xsl:value-of select="@Z"/></xsl:template>
	
	<xsl:template match="WedgeInfo">A1 <xsl:value-of select="@Angle1"/> A2 <xsl:value-of select="@Angle2"/> </xsl:template>
	

</xsl:stylesheet>