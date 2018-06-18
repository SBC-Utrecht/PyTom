<?xml version="1.0" encoding="UTF-8"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

	<xsl:template match="/ParticleList">
		
		<html>
			<title>ParticleList</title>
			<body>
				
				<table border="1">
					<tr bgcolor="#9acd32"><th>Filename</th><th>Rotation (Z1,X,Z2)</th><th>Shift (X,Y,Z)</th><th>Wedge</th></tr>
			
					<xsl:for-each select="Particle">
						<tr>
							<td><xsl:value-of select="@Filename"/></td>
							<td><xsl:apply-templates select="Rotation"/></td>
							<td><xsl:apply-templates select="Shift"/></td>
							<td><xsl:apply-templates select="WedgeInfo"/></td>
						</tr>
					</xsl:for-each>
				</table>
				
			</body>
		</html>
		
	</xsl:template>
	
	<xsl:template match="Rotation">
		<xsl:value-of select="@Z1"/>; 
		<xsl:value-of select="@Z2"/>; 
		<xsl:value-of select="@X"/>
	</xsl:template>
	
	<xsl:template match="Shift">
		<xsl:value-of select="@X"/>; 
		<xsl:value-of select="@Y"/>; 
		<xsl:value-of select="@Z"/>
	</xsl:template>
	
	<xsl:template match="WedgeInfo">
		<xsl:value-of select="@Angle"/>
	</xsl:template>
	

</xsl:stylesheet>