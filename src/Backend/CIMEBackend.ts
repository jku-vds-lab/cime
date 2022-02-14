export class CIMEBackend {
  protected smiles_cache = {};
  protected smiles_highlight_cache = {};
  protected cache = {};

  constructor(
    public readonly baseUrl: string,
    public readonly fetchParams: Parameters<typeof fetch>[1] = {}
  ) {}

  protected handleSmilesCache = (smiles: string, highlight = false) => {
    //already downloaded this image -> saved in smiles cache
    if (highlight) {
      return this.smiles_highlight_cache[smiles];
    } else {
      return this.smiles_cache[smiles];
    }
  };

  protected setSmilesCache = (smiles, highlight = false, data) => {
    if (highlight) this.smiles_highlight_cache[smiles] = data;
    else this.smiles_cache[smiles] = data;
  };

  protected async_cache = async (cached_data) => {
    return cached_data;
  };

  handleCache = (key) => {
    if (this.cache[key]) return Object.assign(this.cache[key]); // return copy of cached object
    return null;
  };

  setCache = (key, value) => {
    this.cache[key] = value;
  };

  handleErrors = (response) => {
    if (!response.ok) {
      throw Error(response.statusText);
    }
    return response;
  };
  handleJSONErrors = (data) => {
    if (Object.keys(data).includes("error")) {
      alert(data["error"]);
    }
    return data;
  };

  public deleteFile = async (filename): Promise<{ deleted: boolean }> => {
    let path = this.baseUrl + "/delete_file/" + filename;

    return fetch(path, {
      ...this.fetchParams,
      method: "GET",
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        alert("file could not be deleted. please, try again");
        console.log(error);
      });
  };

  public getFiles = async (): Promise<{ name: string; id: number }[]> => {
    let path = this.baseUrl + "/get_uploaded_files_list";

    return fetch(path, {
      ...this.fetchParams,
      method: "GET",
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        // alert("could not load uploaded filenames.")
        console.log(error);
      });
  };

  public getDifferenceHighlight = async (
    smilesA: any,
    smilesB: any,
    controller
  ) => {
    const formData = new FormData();
    formData.append("smilesA", smilesA);
    formData.append("smilesB", smilesB);

    let path = this.baseUrl + "/get_difference_highlight";
    let my_fetch;
    if (controller) {
      my_fetch = fetch(path, {
        ...this.fetchParams,
        method: "POST",
        body: formData,
        signal: controller.signal,
      });
    } else {
      my_fetch = fetch(path, {
        ...this.fetchParams,
        method: "POST",
        body: formData,
      });
    }

    return my_fetch
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        console.log(data);
        return data["data"];
      })
      .catch((error) => {
        // alert("could not load structure");
        console.log(error);
      });
  };

  public getStructureFromSmiles = (
    id: string | number,
    smiles: string,
    highlight = false,
    controller
  ) => {
    const cached_data = this.handleSmilesCache(smiles, highlight);
    if (cached_data) {
      return this.async_cache(cached_data);
    }

    const formData = new FormData();
    formData.append("smiles", smiles);
    formData.append("filename", id?.toString());

    let path = this.baseUrl + "/get_mol_img";
    if (highlight) {
      path += "/highlight";
    }

    return fetch(path, {
      ...this.fetchParams,
      method: "POST",
      body: formData,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        this.setSmilesCache(smiles, highlight, data["data"]);
        return data["data"];
      })
      .catch((error) => {
        // alert("could not load structure");
        console.log(error);
      });
  };

  public getStructuresFromSmilesList = async (
    formData: FormData,
    controller?
  ) => {
    return fetch(this.baseUrl + "/get_mol_imgs", {
      ...this.fetchParams,
      method: "POST",
      body: formData,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        if (data["error_smiles"].length > 0) {
          alert(
            "following smiles couldn not be parsed: " + data["error_smiles"]
          );
        }
        return data;
      })
      .catch((error) => {
        if (error.name === "AbortError") {
          console.log("Fetch aborted");
        } else {
          alert("could not load structures");
          console.log(error);
        }
      });
  };

  public getMCSFromSmilesList = async (formData: FormData, controller?) => {
    let my_fetch;
    if (controller) {
      my_fetch = fetch(this.baseUrl + "/get_common_mol_img", {
        method: "POST",
        body: formData,
        signal: controller?.signal,
      });
    } else {
      my_fetch = fetch(this.baseUrl + "/get_common_mol_img", {
        method: "POST",
        body: formData,
      });
    }
    return my_fetch
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((response) => response["data"])
      .catch((error) => {
        // alert("could not get maximum common substructure")
        console.log(error);
      });
  };

  public getSubstructureCount = async (smiles_list, filter) => {
    const formData = new FormData();
    formData.append("smiles_list", smiles_list);
    formData.append("filter_smiles", filter);
    return fetch(this.baseUrl + "/get_substructure_count", {
      method: "POST",
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        if (Object.keys(data).includes("substructure_counts"))
          return data["substructure_counts"];
        else throw Error("Backend responded with error: " + data["error"]);
      })
      .catch((error) => {
        alert("could not find substructure match");
        console.log(error);
      });
  };

  public upload_sdf_file = async (
    file,
    controller?
  ): Promise<{ name: string; id: number }> => {
    // upload the sdf file to the server
    // the response is a unique filename that can be used to make further requests
    const formData_file = new FormData();
    formData_file.append("myFile", file);
    formData_file.append("file_size", file.size);

    const promise = fetch(this.baseUrl + "/upload_sdf", {
      ...this.fetchParams,
      method: "POST",
      body: formData_file,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        if (error.name === "AbortError") {
          console.log("Fetch aborted");
        } else {
          alert("error when uploading file. it might be too big");
          console.log(error);
        }
      });
    return promise;
  };

  public getRepresentationList = async (
    refresh = false,
    id: string | number,
    controller: AbortController
  ) => {
    if (!refresh) {
      const cached_data = this.handleCache("representation_list_" + id);
      if (cached_data && cached_data["rep_list"].length > 0) {
        return this.async_cache(cached_data);
      }
    }
    let path = this.baseUrl + `/get_atom_rep_list/${id}`;

    return fetch(path, {
      ...this.fetchParams,
      method: "GET",
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        this.setCache("representation_list_" + id, data);
        return data;
      })
      .catch((error) => {
        // alert("error when loading representation list")
        console.log(error);
      });
  };

  public calculateHDBScanClusters = async (
    X,
    min_cluster_size,
    min_cluster_samples,
    allow_single_cluster
  ) => {
    const formData = new FormData();
    formData.append("min_cluster_size", min_cluster_size);
    formData.append("min_cluster_samples", min_cluster_samples);
    formData.append("allow_single_cluster", allow_single_cluster);
    formData.append("X", X);
    return fetch(this.baseUrl + "/segmentation", {
      method: "POST",
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        alert("error when calculating clusters");
        console.log(error);
      });
  };
}

// Use the environment variables defined in the .env file
if (!process.env.REACT_APP_CIME_BACKEND_URL) {
  console.error("The ENV-variable REACT_APP_CIME_BACKEND_URL must be set.");
}

export const CIMEBackendFromEnv = new CIMEBackend(
  process.env.REACT_APP_CIME_BACKEND_URL,
  {
    credentials: (process.env.REACT_APP_CIME_BACKEND_CREDENTIALS as RequestCredentials) || 'same-origin',
  }
);
