import React, { useState } from 'react';
import logo from './logo.svg';
import './App.css';
import { tryIt, Application, PSEContextProvider, API } from 'projection-space-explorer';

function App() {
  const [context, setContext] = useState(new API())

  return <PSEContextProvider context={context}>
    <Application />
  </PSEContextProvider>
}

export default App;
