import React, { useState } from 'react';
import logo from './logo.svg';
import './App.css';
import { PSEContextProvider, API, Application, PluginRegistry } from 'projection-space-explorer';
import { ChemPlugin } from './plugins/Cime/ChemPlugin';

PluginRegistry.getInstance().registerPlugin(new ChemPlugin())

function App() {
  const [context, setContext] = useState(new API())

  return <PSEContextProvider context={context}>
    <Application config={{}} />
  </PSEContextProvider>
}

export default App;
