import { useState } from 'react';
import './App.css';
import { PSEContextProvider, API, Application, PluginRegistry } from 'projection-space-explorer';
import { ChemPlugin } from './plugins/Cime/ChemPlugin';
import { DatasetTabPanel } from './DatasetTabPanel';
import { CimeAppBar } from './CimeAppBar';


export const DEMO = false

/**
 *               <a href={"https://jku-vds-lab.at"} target={"_blank"}><VDSLogo style={{ height: 48, width: 48 }}></VDSLogo></a>
              {frontend_utils.CHEM_PROJECT && <a href={"https://www.bayer.com"} target={"_blank"}><img style={{ height: 48, marginLeft: 48 }} src={"textures/bayer-logo.svg"} alt="Powered By Bayer" /></a>}
              <Typography variant="h6" style={{ marginLeft: 48, color: "rgba(0, 0, 0, 0.54)" }}>
                {frontend_utils.CHEM_PROJECT ? "CIME: Chem-Informatics Model Explorer" : "Projection Space Explorer"}
              </Typography>
 */

PluginRegistry.getInstance().registerPlugin(new ChemPlugin())


function App() {
  const [context, setContext] = useState(new API())

  return <PSEContextProvider context={context}>
    <Application
      config={{
        preselect: { url: 'datasets/test.sdf' }
      }}
      features={{
        disableEmbeddings: {
          tsne: true,
          forceatlas: true
        }
      }}
      overrideComponents={{
        datasetTab: DatasetTabPanel,
        appBar: CimeAppBar
      }} />
  </PSEContextProvider>
}

export default App;