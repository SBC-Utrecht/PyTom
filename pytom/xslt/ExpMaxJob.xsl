<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
            
            <xsl:template match="/ExpectationMaximisationJob">
                
                <html>
					<style type="text/css"> 
					 h1[type="Main"] { color:#66CCFF; font-size:xx-large; font-family:"Verdana" }
					 h1[type="SubMain"] { color:white; font-size:large; font-family:"Verdana" }  
					 h5[type="Copyright"] { color:white; font-size:x-small; font-family:Times vertical-align: bottom;} 
					 p[type="MainText"]{font-size:medium; color:white; font-family:"Verdana"} 
					 TR{font-family:"Verdana"; color:#FFFFFF;}
					 TH{font-family:"Verdana"; color:#FFFFFF;}
					 TD[type="TableCaption"]{font-family:"Verdana"; color:#FFFFFF; font-weight:bold;}
					 TD[type="TableRest"]{font-family:"Verdana"; color:#FFFFFF;}
					 a { color:#33FF33; font-size:medium; font-family:"Verdana" }
					 body { background-color:black; }
					 *[lang=en]{quotes:"'" "'";} 
					</style>
                    <title>Expectation Maximisation Job</title>
                    <body>
                        <center>
                        <h1 type="Main">Expectation Maximisation Job</h1>
                        <table border="1">
                            <xsl:apply-templates select="Description"/>
                            <tr><td bgcolor="#AAAADD" TITLE="Current alignment reference">Start reference</td><td><xsl:apply-templates select="Description/Reference"/></td></tr>
                            <tr><td bgcolor="#AAAAA0" TITLE="Alignemnt mask">Mask</td><td><xsl:apply-templates select="Description/Mask"/></td></tr>
                            <tr><td bgcolor="#AAAADD" TITLE="Method used for scoring">Score</td><td><xsl:apply-templates select="Description/Score"/></td></tr>
                            <tr><td bgcolor="#AAAAA0">Angles</td><xsl:apply-templates select="Description/Angles"/></tr>
                            <tr><td bgcolor="#AAAADD" TITLE="Detailed information about your sample in Angstrom">Sample Information</td><td><xsl:apply-templates select="Description/SampleInformation"/></td></tr>
                            <tr><td bgcolor="#AAAAA0" TITLE="Symmetry. 1 means no symmetry, 2 is two fold symmetry around symmetry axis. Z is the default">Symmetry</td><td><xsl:apply-templates select="Description/Symmetry"/></td></tr>
							<xsl:apply-templates select="Description/Preprocessing"/>    
                        </table>                
                        <xsl:apply-templates select="Description/ParticleList"/>
                        </center>
                    </body>
                    
                </html>
            </xsl:template>
            
            <xsl:template match="Description">
                <tr><td bgcolor="#AAAADD" TITLE="Where are the results stored?">Destination</td><td><xsl:value-of select="@Destination"/></td></tr>
                <tr><td bgcolor="#AAAAA0" TITLE="How many iterations">Number Iterations</td><td><xsl:value-of select="@NumberIterations"/></td></tr>
                <tr><td bgcolor="#AAAADD" TITLE="">Include only n* numberParticles with the highest scores for averaging.</td><td><xsl:value-of select="@ResultClassification"/>  </td></tr>
                <tr><td bgcolor="#AAAAA0" TITLE="How often is the data shrunken? libtomc notation!">Binning</td><td><xsl:value-of select="@Binning"/></td></tr>
                <tr><td bgcolor="#AAAADD" TITLE="FSCCriterion. Specify at which FSC the resolution will be determined.">FSCCriterion</td><td><xsl:value-of select="@FSCCriterion"/></td></tr>
                <tr><td bgcolor="#AAAAA0" TITLE="AdaptiveResolution. Enable / disable adaptive lowpass and refinement angle to current resolution.">AdaptiveResolution</td><td><xsl:value-of select="@AdaptiveResolution"/></td></tr>
            </xsl:template>
            
            <xsl:template match="Reference">
                <xsl:value-of select="@File"/><br/>
            </xsl:template>
            
            <xsl:template match="Mask">
                <xsl:value-of select="@Filename"/><br/>
            </xsl:template>
            
            <xsl:template match="Score">
                <xsl:value-of select="@Type"/><br/>
            </xsl:template>
            
            <xsl:template match="Angles">
                <xsl:choose>
                    <xsl:when test="@Type='FromEMFile'"><td TITLE="Angle list stored on disk">FromEMFile : <xsl:value-of select="@File"/></td></xsl:when>
                    <xsl:when test="@Type='Equidistant'"><td TITLE="Angle list generated around a predifined rotation. See F. Foerster PhD. Thesis">Equidistant : Increment <xsl:value-of select="@Increment"/>  Shells<xsl:value-of select="@Shells"/></td></xsl:when>
                </xsl:choose>
            </xsl:template>
            
            <xsl:template match="Parameters">
                Increment <xsl:value-of select="Increment"/>,
                Shells <xsl:value-of select="Shells"/>
            </xsl:template>
            
            <xsl:template match="SampleInformation">
                PixelSize <xsl:value-of select="@PixelSize"/>, Particle Diameter <xsl:value-of select="@ParticleDiameter"/><br/>
            </xsl:template>
            
            <xsl:template match="Symmetry">
                <xsl:value-of select="@NumberSymmetries"/>
                Symmetry Axis Rotation: X <xsl:value-of select="@Theta"/>, Z2 <xsl:value-of select="@Psi"/> (0,0 defaults to Z as symmetry axis)<br/>
            </xsl:template>
            
            <xsl:template match="Preprocessing">
                 <xsl:apply-templates select="Bandpass"/>
            </xsl:template>
            
            <xsl:template match="Bandpass">
                 <tr><td bgcolor="#AAAADD" TITLE="Sliding bandpass set to lowest / highest frequency (in Nyquist) and smoothing (in pixels)">Bandpass</td><td>Lowest F <xsl:value-of select="@LowestFrequency"/>, Highest F <xsl:value-of select="@HighestFrequency"/>, Smooth <xsl:value-of select="@Smooth"/></td></tr>
            </xsl:template>
            
            <xsl:template match="ParticleList">
                 <br/><br/>
                <table border="1">
                    <tr bgcolor="#9acd32"><th TITLE="Absolute path to particle ">Filename</th><th TITLE="Rotation in ZXZ paradigm. Z1 : first rotation around Z , X rotation around X , Z2 second rotation around Z ">Rotation (Z1,X,Z2)</th><th TITLE="Shift of particle into the center">Shift (X,Y,Z)</th><th TITLE="Wedge size of particle">Wedge</th></tr>
                
                    <xsl:for-each select="Particle">
                        <tr>
                            <td><xsl:value-of select="@Filename"/></td>
                            <td><xsl:apply-templates select="Rotation"/></td>
                            <td><xsl:apply-templates select="Shift"/></td>
                            <td><xsl:apply-templates select="WedgeInfo"/></td>
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