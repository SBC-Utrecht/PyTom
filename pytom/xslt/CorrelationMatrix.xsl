<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
            
            <xsl:template match="/CorrelationMatrixJob">
                <html>
					<style type="text/css"> 
						 h1[type="Main"] { color:#66CCFF; font-size:xx-large; font-family:"Verdana" }
						 h1[type="SubMain"] { color:white; font-size:large; font-family:"Verdana" }  
						 h5[type="Copyright"] { color:white; font-size:x-small; font-family:Times vertical-align: bottom;}  
						 TR{font-family:"Verdana"; color:#FFFFFF;}
						 TH{font-family:"Verdana"; color:#FFFFFF;}
						 TD[type="TableCaption"]{font-family:"Verdana"; color:#FFFFFF; font-weight:bold;}
						 TD[type="TableRest"]{font-family:"Verdana"; color:#FFFFFF;}
						 a { color:#33FF33; font-size:medium; font-family:"Verdana" }
						 body { background-color:black; }
						 *[lang=en]{quotes:"'" "'";} 
					</style>
                    <title>Correlation Matrix Job</title>
                    <body>
                        <center>
                        <h1 type="Main">Correlation Matrix Job</h1>
                        <table border="1">
                        <tr><td bgcolor="#AAAADD" TITLE="File where matrix will be stored" type="TableCaption">Correlation Matrix Name</td><td type="TableRest"><xsl:value-of select="@ResultMatrixName"/></td></tr>
                        <tr><td bgcolor="#AAAADD" TITLE="Will wedge be applied?" type="TableCaption">Apply Wedge</td><td type="TableRest"><xsl:value-of select="@ApplyWedge"/></td></tr>
                        <tr><td bgcolor="#AAAADD" TITLE="Lowest and highest frequency" type="TableCaption">Bandpass</td><td type="TableRest"><xsl:value-of select="@LowestFrequency"/>;<xsl:value-of select="@HighestFrequency"/></td></tr>
                        <tr><td bgcolor="#AAAADD" TITLE="How many pixels are combined for downsampling. (1 means downsampling disabled)" type="TableCaption">Downsampling size</td><td type="TableRest"><xsl:value-of select="@Binning"/></td></tr>
                        </table>
                        <xsl:apply-templates select="ParticleList"/>
                        </center>
                    </body>
                </html>
                
                
            </xsl:template>
            
            <xsl:template match="ParticleList">
                 <br/><br/>
                <table border="1">
                    <tr bgcolor="#9acd32"><th TITLE="Absolute path to particle ">Filename</th><th TITLE="Rotation in ZXZ paradigm. Z1 : first rotation around Z , X rotation around X , Z2 second rotation around Z ">Rotation (Z1,X,Z2)</th><th TITLE="Shift of particle into the center">Shift (X,Y,Z)</th><th TITLE="Wedge size of particle">Wedge</th></tr>
                
                    <xsl:for-each select="Particle">
                        <tr>
                            <td type="TableRest"><xsl:value-of select="@Filename"/></td>
                            <td type="TableRest"><xsl:apply-templates select="Rotation"/></td>
                            <td type="TableRest"><xsl:apply-templates select="Shift"/></td>
                            <td type="TableRest"><xsl:apply-templates select="WedgeInfo"/></td>
                        </tr>
                    </xsl:for-each>
                </table>
            </xsl:template>
            
            <xsl:template match="Rotation">
                <xsl:value-of select="@Z1"/>;
                <xsl:value-of select="@X"/>;
                <xsl:value-of select="@Z2"/>
            </xsl:template>
            
            
            <xsl:template match="Shift">
                <xsl:value-of select="@X"/>; 
                <xsl:value-of select="@Y"/>; 
                <xsl:value-of select="@Z"/>
            </xsl:template>
            
</xsl:stylesheet>