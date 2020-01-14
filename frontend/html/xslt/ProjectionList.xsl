<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:template match="/ProjectionList">
		<table class="projectionList">
					<xsl:for-each select="Projection">
						
						<tr class="projectionList">
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:value-of select="@Filename"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							<xsl:element name="b">
								Pr<xsl:value-of select="position()"/>
							</xsl:element> 
						</xsl:element>
							
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:apply-templates select="@TiltAngle"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							Tilt
						</xsl:element>			
							
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:apply-templates select="@OffsetX"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							OffsetX
						</xsl:element>	
							
						<xsl:element name="td">
							<xsl:attribute name="class">innerList</xsl:attribute>
							<xsl:attribute name="align">center</xsl:attribute>
							<xsl:attribute name="onMouseOver">document.getElementById('Help').innerHTML='<xsl:apply-templates select="@OffsetY"/>';</xsl:attribute>
							<xsl:attribute name="onMouseOut">document.getElementById('Help').innerHTML='';</xsl:attribute>
							OffsetY
						</xsl:element>
							
						</tr>
					</xsl:for-each>
				</table>
	</xsl:template>
</xsl:stylesheet>