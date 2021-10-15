import { useState } from 'react';
import { PSEContextProvider, API, Application, PluginRegistry, RootState, createRootReducer } from 'projection-space-explorer';
import { ChemPlugin } from './Cime/ChemPlugin';
import { DatasetTabPanel } from './DatasetTabPanel';
import { CimeAppBar } from './CimeAppBar';
import { LineUpContext } from './LineUpContext';
import {setLineUpInput_visibility} from 'projection-space-explorer';


export const DEMO = false

console.log(setLineUpInput_visibility)

const reducers = {
  lineUpInpu: null
}


PluginRegistry.getInstance().registerPlugin(new ChemPlugin())

//PluginRegistry.getInstance().registerReducer(lineUpInput)
//PluginRegistry.getInstance().registerReducer()
export type AppState = RootState & {

}


function App() {
  const [context, setContext] = useState(new API<AppState>(null, createRootReducer(reducers)))

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