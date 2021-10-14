import { useState } from 'react';
import { PSEContextProvider, API, Application, PluginRegistry } from 'projection-space-explorer';
import { ChemPlugin } from './Cime/ChemPlugin';
import { DatasetTabPanel } from './DatasetTabPanel';
import { CimeAppBar } from './CimeAppBar';
import { LineUpContext } from './LineUpContext';
import lineUpInput from 'projection-space-explorer';


export const DEMO = false

PluginRegistry.getInstance().registerPlugin(new ChemPlugin())

//PluginRegistry.getInstance().registerReducer(lineUpInput)
//PluginRegistry.getInstance().registerReducer()


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
        appBar: CimeAppBar,
        detailViews: [{
          name: 'lineup',
          //@ts-ignore
          view: LineUpContext
        }]
      }} />
  </PSEContextProvider>
}

export default App;